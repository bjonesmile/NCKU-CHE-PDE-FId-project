from gams import *
import sys
import os
import threading
import numpy as np
import itertools
import math

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import csv
import dill as pickle

class Temp_distribute1D():
    def __init__(self,dt,dx):
        self.total_sim_time = 350
        self.x_grid_num = int(1/dx)

        self.ws = GamsWorkspace(working_directory = os.getcwd())
        self.dt = dt
        self.dx = dx
        self.BC_L = 1.0 #B.C. x = x1
        self.BC_R = 1.0 #B.C. x = x20
        self.IC_X = np.zeros((self.x_grid_num)) #I.C. t = t1
        self.temp_lb = 0.25
        self.drawing = 1
        
        self.plot_t = 1
        self.look_point = np.arange(1,21,1,dtype=np.int)

        self.jobGAMSfile = "temp_slab_FId.gms"

    def check_region_x(self,set_x):
        if int(set_x[1:]) in self.look_point:
            return True
        return False

    def check_region_t(self,set_t):
        if int(set_t[1:]) % self.plot_t == 0 :
            return True
        return False

    def solve_result(self):
        db = self.ws.add_database()
        job = self.ws.add_job_from_file(self.jobGAMSfile)
        opt = self.ws.add_options()
        opt.defines["gdxincname"] = db.name
        opt.all_model_types = "sbb"
        opt.idcgdxoutput = "test.gdx"
        dt = db.add_parameter("dt", 0)
        dt.add_record().value = self.dt
        dx = db.add_parameter("dx", 0)
        dx.add_record().value = self.dx
        
        # const linear phasor step

        init_state = db.add_parameter("init_state", 1, "I.C. for each x grid point in time horizon for uncertain parameter temperature")        
        self.set_IC_X(0.35,1,3,ctype='const')
        for i in range(1,int(1/self.dx)+1):
            init_state.add_record('x'+str(i)).value = self.IC_X[i-1]

        bc_l = db.add_parameter("BC_L",1,"B.C. for each x1 grid point")
        self.set_BC_L(0.5,1,1,ctype='const')
        for i in range(1,self.total_sim_time+1):
            bc_l.add_record('t'+str(i)).value = self.BC_L[i-1]

        bc_r = db.add_parameter("BC_R",1,"B.C. for each x1 grid point")
        self.set_BC_R(0.5,1,1,ctype='const')
        for i in range(1,self.total_sim_time+1):
            bc_r.add_record('t'+str(i)).value = self.BC_R[i-1]
        
        job.run(opt,databases = db)
        obj = job.out_db.get_variable("obj").first_record().level

        t_len = job.out_db["T"].get_number_records()

        fig, axs = plt.subplots(2,1,sharex=True,sharey=True)
 
        x_grid = np.linspace(0,1,self.x_grid_num)
        t_grid = np.linspace(0,self.total_sim_time,self.total_sim_time)
        data = np.empty((t_len//self.plot_t,len(self.look_point)))
        print(data.shape)
        i, j = 0, 0 
        for rec in job.out_db["u"] :
            if self.check_region_t(rec.key(1)) and self.check_region_x(rec.key(0)):
                data[i][j] = rec.level
                i += 1
                if i >= data.shape[0] :
                    i = 0
                    j += 1
        display_time_list = np.arange(30,350,40)
        for t in display_time_list :
            axs[0].plot(x_grid,data[t][:],label=("t="+str(round((t)*self.plot_t*self.dt,3))))

        major_ticks = np.arange(0,101,10)/100
        #minor_ticks = np.arange(0,101,20)/100
        axs[0].set_xticks(major_ticks)
        #ax.set_xticks(minor_ticks, minor=True)
        axs[0].grid(which='both',color='y', linestyle='--', linewidth=1, alpha=0.3)
        axs[0].legend(loc='best')

        display_time_list = np.arange(1,11,dtype=np.int)
        for t in display_time_list :
            axs[1].plot(x_grid,data[t][:],label=("t="+str(round((t)*self.plot_t*self.dt,3))))
        axs[1].set_xticks(major_ticks)
        axs[1].grid(which='both',color='y', linestyle='--', linewidth=1, alpha=0.3)
        axs[1].legend(loc='best')
        fig.tight_layout()
        plt.show()
        plt.close()

        X, Y = np.meshgrid(x_grid,t_grid)

        def plotfunc(plotdata):
            axis = plt.gca()
            surface = axis.plot_surface(plotdata['xs'],plotdata['ys'],plotdata['zs'], rstride=1, cstride=1, cmap='coolwarm')            
            if len(plt.gcf().axes) == 1 :
                plt.colorbar(surface,shrink=1.0,aspect=20)
            axis.set_title("temperature distribute", fontsize=16)
            axis.set_xlabel("X grid point",fontsize=16)
            axis.set_ylabel("T grid point",fontsize=16)
            

        fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
        plotdata = {'xs': X, 'ys': Y, 'zs': data, 'plotfun': plotfunc}
        plotdata['plotfun'](plotdata)
        with open("test.pickle",'wb') as file :
            pickle.dump(plotdata, file)
        plt.show()
        plt.close()
        return obj

    def load_fig(self, picklefile):
        plt.figure().add_subplot(111, projection='3d')
        plotdata = pickle.load(open(picklefile, 'rb'))
        plotdata['plotfun'](plotdata)
        plt.show()

    def set_IC_X(self,magnitude,theta,omega,ctype='const'):
        if ctype == 'const' :
            ary = np.ones(self.total_sim_time)
            ary = ary*magnitude
        elif ctype == 'linear' :
            ary = np.linspace(0.5,magnitude,350)
        elif ctype == 'phasor':
            ary = np.linspace(0,omega*math.pi,num=self.x_grid_num,endpoint=True)
            ic = magnitude*np.cos(ary+theta*math.pi/2)
            ic = abs(ic)+self.temp_lb
            ary = ic
        
        self.IC_X = ary

    def set_BC_R(self,magnitude,theta,omega,ctype='const'):
        if ctype == 'const' :
            ary = np.ones(350)
            ary = ary*magnitude
        elif ctype == 'linear' :
            ary = np.linspace(0.5,magnitude,350)
        elif ctype == 'phasor':
            ary = np.linspace(0,omega*math.pi,num=self.total_sim_time,endpoint=True)
            bc = magnitude*np.cos(ary+theta*math.pi/2)
            bc = abs(bc)
            ary = bc
        elif ctype == 'step':
            step_num = 10
            ary = np.ones(self.total_sim_time)
            step_set = self.total_sim_time//step_num
            for i in range(step_num):
                ary[i*step_set:(i+1)*step_set] = magnitude+0.1*i
                if i == step_num-1:
                    ary[i*step_set:] = magnitude+0.1*i
            t = np.arange(0,350)
            plt.plot(t,ary)
            plt.show()
        self.BC_R = ary

    def set_BC_L(self,magnitude,theta,omega,ctype='const'):
        if ctype == 'const' :
            ary = np.ones(350)
            ary = ary*magnitude
        elif ctype == 'linear' :
            ary = np.linspace(0.5,magnitude,350)
        elif ctype == 'phasor':
            ary = np.linspace(0,omega*math.pi,num=self.total_sim_time,endpoint=True)
            bc = magnitude*np.cos(ary+theta*math.pi/2)
            bc = abs(bc)
            ary = bc
        elif ctype == 'step':
            step_num = 10
            ary = np.ones(self.total_sim_time)
            step_set = self.total_sim_time//step_num
            for i in range(step_num):
                ary[i*step_set:(i+1)*step_set] = magnitude+0.1*i
                if i == step_num-1:
                    ary[i*step_set:] = magnitude+0.1*i
            t = np.arange(0,350)
            plt.plot(t,ary)
            plt.show()
        self.BC_L = ary

    def set_temp_base_value(self,temp):
        self.temp_lb = temp
        return

if __name__ == '__main__' :
    test = Temp_distribute1D(0.001,0.05)
    test.set_temp_base_value(0.25)
    result = test.solve_result()
    test.load_fig('test.pickle')

    print(result)

    exit(0)