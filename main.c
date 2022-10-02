#include "dynamics.h"
#include <chrono>
// #include "udp.h"
#include <string>
#include <thread>
#include <signal.h>
#include <sys/time.h>
#include "mytimer.h"
#include <signal.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

using namespace std;
using namespace RigidBodyDynamics;
using namespace RigidBodyDynamics::Math;

const char *modelfile = "./webserver/model/effort/ER210-stl.urdf";
Model *model = NULL;
int dofs = 6;
VectorNd q0, q1, q, qik, qdot, qddot, tau, tauid, q_last, qdot_last, qtemp;
Vector3d endpos0, endrpy0, endpos1, endrpy1, endpos, endrpy;
Matrix3d endori,endori0,endori1,ori_ik;
Vector3d pos_ik;
ctrl pid;

double tt = 2.0, dt = 0.001, timeused = 0;
double errkine = 0;
double errdyna = 0;
double amp = 0.1, w = 1, t = 0;
VectorNd errMax;
dynamics rbdl;
// timer
size_t timer = 0;
double starttime, starttime_1;
chrono::system_clock::time_point tm_start;
bool toQuit = false;
int test_mode = 0;
bool remote = true;
// void sendinfo(const char head[4], double data[], int len);
bool tosend = false;


//***********************wsserver**************************************
#include <ws.h>
void ws_sendinfo(double data[], int len)
{
    float fdata[len];
    for (int i = 0; i < len; i++)
        fdata[i] = data[i];
    ws_sendframe_bin(NULL, (char *)fdata, len*4);
}


void ws_onopen(ws_cli_conn_t *client)
{
	char *cli;
	cli = ws_getaddress(client);
#ifndef DISABLE_VERBOSE
	printf("Connection opened, addr: %s\n", cli);
#endif
}

/**
 * @brief Called when a client disconnects to the server.
 *
 * @param client Client connection. The @p client parameter is used
 * in order to send messages and retrieve informations about the
 * client.
 */
void ws_onclose(ws_cli_conn_t *client)
{
	char *cli;
	cli = ws_getaddress(client);
#ifndef DISABLE_VERBOSE
	printf("Connection closed, addr: %s\n", cli);
#endif
}


void ws_onmessage(ws_cli_conn_t *client,
	const unsigned char *msg, uint64_t size, int type)
{
	char *cli;
	cli = ws_getaddress(client);
#ifndef DISABLE_VERBOSE
	// printf("I receive a message: %s (size: %" PRId64 ", type: %d), from: %s\n",
	// 	msg, size, type, cli);
#endif

	/**
	 * Mimicks the same frame type received and re-send it again
	 *
	 * Please note that we could just use a ws_sendframe_txt()
	 * or ws_sendframe_bin() here, but we're just being safe
	 * and re-sending the very same frame type and content
	 * again.
	 *
	 * Client equals to NULL: broadcast
	 */
    if(strncmp((const char*)msg, "data", 4) == 0)
    {
        // tosend = true;
        ws_sendinfo((double *)rbdl.q.data(), 6);
    }
	//
}

/**
 * @brief Main routine.
 *
 * @note After invoking @ref ws_socket, this routine never returns,
 * unless if invoked from a different thread.
 */
int ws_init(void)
{
	struct ws_events evs;
	evs.onopen    = &ws_onopen;
	evs.onclose   = &ws_onclose;
	evs.onmessage = &ws_onmessage;
	ws_socket(&evs, 8080, 1, 1000); /* Never returns. */

	/*
	 * If you want to execute code past ws_socket, invoke it like:
	 *   ws_socket(&evs, 8080, 1, 1000)
	 */

	return (0);
}
//***********************wsserver**************************************

double rsttm(void)
{
    tm_start = chrono::system_clock::now();
    return 0;
}
double gettm(void)
{
    chrono::duration<double> elapsed = chrono::system_clock::now() - tm_start;
    return elapsed.count();
}

// void sendinfo(const char head[4], float data[], int len)
// {
//     char buf[len * 4 + 4];
//     sprintf(buf, head, 4);
//     for (int i = 0; i < len; i++)
//     {
//         memcpy(&buf[4] + i * 4, &data[i], 4);
//     }
//     udp_send(buf, len * 4 + 4);
// }

// void sendinfo(const char head[4], double data[], int len)
// {
//     float fdata[len];
//     for (int i = 0; i < len; i++)
//         fdata[i] = data[i];
//     sendinfo(head, fdata, len);
// }



// void on_remoteCmd(char *cmd)
// {
//     // execCmd(cmd[0]);
//     printf("remote cmd:%s\n", cmd);
// }

struct state_machine
{
    int episoid_step_count = 1000;
    int episoid_step_idx = 0;
    int episoid_idx = 0;
    bool stopped = false;
} SM;

void stopTest()
{
    if(test_mode==0)stop_timer(timer);
    SM.episoid_idx = 0;
    SM.episoid_step_idx = 0;
}

void reset_ctrl()
{
    pid.Kp = VectorNd::Ones(dofs);
    pid.Kp << 50000.0,50000.0,25000.0,100.0,100.0,100.0;
    pid.Ki = 0.0*VectorNd::Ones(dofs);
    pid.Kd = 2.0*VectorNd::Ones(dofs);
    pid.Kd << 1500,1500,1500,2,2,2;
    pid.err = VectorNd::Zero(dofs);
    pid.err_d = VectorNd::Zero(dofs);
    pid.err_i = VectorNd::Zero(dofs);
    pid.u = VectorNd::Zero(dofs);
    pid.err_last =  VectorNd::Zero(dofs);
}

void SM_proc(size_t timer_id, void *user_data)
{
	static int batch_count = 0;
    static int tick_count = 0;

    // printf("state machine episoid:%d,step:%d\n",SM.episoid_idx,SM.episoid_step_idx);
    switch (SM.episoid_idx)
    {
    case 0: // FK test
        if (SM.episoid_step_idx == 0)
        {
            q1[0] = 0 * M_PI / 4;
            q1[1] = M_PI*1 / 4;
            q1[2] = 0*M_PI / 2;
            q1[4] = 1*M_PI / 4;
            SM.episoid_step_count = 1000;
            starttime = gettm();
            starttime_1 = gettm();
            if(test_mode != 2)ws_sendframe_txt(NULL, "test");
            //rbdl.printLinksCom();
        }
        else
        {
            q = q0 + (q1 - q0) * SM.episoid_step_idx / SM.episoid_step_count;
            rbdl.q = q;
            rbdl.forwardKinematics(q, endpos, endrpy);
        }
        SM.episoid_step_idx++;
        if (SM.episoid_step_idx == SM.episoid_step_count)
        {
            SM.episoid_idx ++;
            SM.episoid_step_idx = 0;
            timeused = gettm() - starttime;
            printf("gohome n=%d, time=%.6f s, single time=%.6f us \n", SM.episoid_step_count, timeused, 1000000 * timeused / SM.episoid_step_count);
        }
        break;
    case 1: // FK test
        if (SM.episoid_step_idx == 0)
        {
            SM.episoid_step_count = 6 * 1000;
            starttime = gettm();
        }
        else
        {
            amp = M_PI / 4;
            w = 1.0f;
            t = SM.episoid_step_idx * dt;
            int ax = SM.episoid_step_idx / 1000;
            q[ax] = q1[ax] + amp * sin(w * 2 * M_PI * t);

            rbdl.forwardKinematics(q, endpos, endori);
            rbdl.q = q;
            // if(remote)sendinfo("qpos", (double *)q.data(), 6);
        }
        SM.episoid_step_idx++;
        if (SM.episoid_step_idx == SM.episoid_step_count)
        {
            SM.episoid_idx++;
            SM.episoid_step_idx = 0;
            timeused = gettm() - starttime;
            printf("FK n=%d, time=%.6f s, single time=%.6f us \n", SM.episoid_step_count, timeused, 1000000 * timeused / SM.episoid_step_count);
        }
        break;
    case 2: // IK test
        if (SM.episoid_step_idx == 0)
        {
            rbdl.setJointPositions(q1);
            rbdl.getEndPose(endpos1, endori1);
            SM.episoid_step_count = 6 * 2000;
            starttime = gettm();
        }
        else
        {
            amp = 0.5f;
            w = 0.5f;
            t = SM.episoid_step_idx * dt;
            int ax = SM.episoid_step_idx / 2000;
            if (ax < 3)
            {
                endpos[ax] = endpos1[ax] + amp * sin(w * 2 * M_PI * t);
            }
            else
            {
                Matrix3d rotmat;
                if (ax == 3)
                    rotmat = rotx(0.5 * M_PI * sin(w * 2 * M_PI * t));
                if (ax == 4)
                    rotmat = roty(0.5 * M_PI * sin(w * 2 * M_PI * t));
                if (ax == 5)
                    rotmat = rotz(0.5 * M_PI * sin(w * 2 * M_PI * t));
                endori = rotmat * endori1;
            }
            rbdl.inverseKinematics(endpos, endori, q, qik);
            rbdl.q = q;
            // if(remote)sendinfo("qpos", (double *)qik.data(), 6);
            q = qik;
        }
        SM.episoid_step_idx++;
        if (SM.episoid_step_idx == SM.episoid_step_count)
        {
            SM.episoid_idx++;
            SM.episoid_step_idx = 0;
            timeused = gettm() - starttime;
            printf("IK n=%d, time=%.6f s, single time=%.6f us \n", SM.episoid_step_count, timeused, 1000000 * timeused / SM.episoid_step_count);
        }
        break;
    case 3://step
        if (SM.episoid_step_idx == 0)
        {
            rbdl.setJointPositions(q1);
            rbdl.getEndPose(endpos1, endori1);
            SM.episoid_step_count = 2000;
            starttime = gettm();
			rbdl.resetDynamics();
            //rbdl.printLinksCom();
            reset_ctrl();
            errMax = pid.err;
            rbdl.qd[0] = 0.1*M_PI;
            rbdl.qd[1] = 0.1*M_PI;
        }
        else
        {
            amp = M_PI / 4;
            w = 1.0f * M_PI;
            t = SM.episoid_step_idx * dt;
            //int i = 1;
            q = q1;

            for(int i=0;i<0;i++)
            {
                q[i] = q1[i] + amp * (1 - cos(w * t));
                qdot[i] = w * amp * sin(w * t) ;
                qddot[i] = w * w  * amp * cos(w * t);
            }

            //rbdl.forwardKinematics(q, endpos, endori);
            pid.err = (q - rbdl.q); 
            pid.err_d = (qdot - rbdl.qd); 
            pid.u = pid.Kp.array()  * pid.err.array()  + pid.Kd.array()  * pid.err_d.array() ;
            pid.err_last =  pid.err;
            //rbdl.inverseDynamics(rbdl.q, rbdl.qd, rbdl.qdd, tauid);
            rbdl.nonlinearEffects(rbdl.q, rbdl.qd, tauid);
            //rbdl.nonlinearEffects(rbdl.q, VectorNd::Zero(dofs), tauid);
            rbdl.tau = 1.0f*tauid + 0.0f*pid.u;
            rbdl.step(dt);
            //rbdl.getJointPositions(q);
            errMax = pid.err.cwiseMax(errMax);
            if(SM.episoid_step_idx % 100==0)
            {
                //cout << "tauid:" << tauid.transpose() << "pid.u:" << pid.u.transpose() << endl;
                // cout << "pid.err:" << rbdl.qdd.transpose() << endl;
            }

        }
        SM.episoid_step_idx++;
        if (SM.episoid_step_idx == SM.episoid_step_count)
        {
            SM.episoid_idx++;
            SM.episoid_step_idx = 0;
            timeused = gettm() - starttime;
            printf("dynaSim n=%d, time=%.6f s, single time=%.6f us \n", SM.episoid_step_count, timeused, 1000000 * timeused / SM.episoid_step_count);
        }

        break;
    case 4:
		batch_count++;
		if(batch_count >= 1)
        {
             SM.stopped =  true;
        }
		else
        {
			SM.episoid_idx = 0;
            SM.episoid_step_idx = 0;
        }
        
        break;
    default:
        break;
    }

    if(SM.stopped ==  true)
    {
        stopTest();
        if(test_mode != 2)ws_sendframe_txt(NULL, "stop");
    }

    // if (remote && (test_mode != 2))
    // {
    //     double td = (gettm() - starttime_1);
        
    //     if ( td > (1.0/60.0))
    //     {
    //         starttime_1 = gettm();
    //         ws_sendinfo((double *)rbdl.q.data(), 6);
    //     }
    // }

    // if(tosend)
    // {
    //     ws_sendinfo((double *)rbdl.q.data(), 6);
    //     tosend = false;
    // }
    tick_count++;

}

void benchtest(int mode)
{

    SM.episoid_idx = 0;
    SM.episoid_step_idx = 0;
    starttime = gettm();

    test_mode = mode;
    SM.stopped =  false;
    if(mode == 0)
    {
        timer = start_timer(1, SM_proc, TIMER_PERIODIC, NULL);
    }
    else
    {
        while(!SM.stopped)
            SM_proc(0,NULL);
    }
}

bool execCmd(char cmd)
{
    bool ret = true;
    switch (cmd)
    {
    case 'Q':
    case 'q':
        toQuit = true;
        break;
    case 's':
        SM.stopped = true;
        break;
    case 'r':
		rbdl.setJointPositions(q0);
		ws_sendinfo((double *)q0.data(), 6);
        break;
    case 'z':
        //     arm.setZeroPos();
        break;
    case 't':
        // demo();
        benchtest(0);
        break;
    case '1':
        // demo();
        benchtest(1);
        break;
    case '2':
        // demo();
        benchtest(2);
        break;
    default:
        cout << "unkown command\n"
             << endl;
        break;
    }
    return ret;
}

void handleInput()
{
    char cmd;
    while (!toQuit)
    {
        char cmdinfo[1024] = " t:\t begin test\n s:\t stop test\n r:\t reset pose\n 1:\t full speed test-I\n 1:\t full speed-II\n q:\t quit\n:";
        cout << cmdinfo;
        cin >> cmd;
        if (!execCmd(cmd))
            break;
        // usleep(200000);
    }
}

void init()
{
    rbdl.loadmodelfile(modelfile);
    endori = Matrix3d::Zero();
    endori0 = Matrix3d::Zero();
    endori1 = Matrix3d::Zero();
    pos_ik = Vector3d::Zero();
    ori_ik = Matrix3d::Zero();
    dofs = rbdl.getDofs();
    q0 = VectorNd::Zero(dofs);
    q1 = VectorNd::Zero(dofs);
    q = VectorNd::Zero(dofs);
    qik = VectorNd::Zero(dofs);
    qdot = VectorNd::Zero(dofs);
    qddot = VectorNd::Zero(dofs);
    tau = VectorNd::Zero(dofs);
    tauid = VectorNd::Zero(dofs);
    endori = Matrix3d::Zero();
    endori0 = Matrix3d::Zero();
    endori1 = Matrix3d::Zero();
    qtemp = VectorNd::Zero(dofs);
	init_timer();
    ws_init();
    system("python3 -m http.server  --directory ./webserver/ 8001 &");
}

void finish()
{
    system("sudo killall python3");
}

int main()
{
    init();
    handleInput();
    finish();
    return 0;
}
