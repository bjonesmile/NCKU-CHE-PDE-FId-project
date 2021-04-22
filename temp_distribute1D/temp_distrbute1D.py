from gams import *
import sys
import os
import threading
import numpy as np
import itertools
import matplotlib.pyplot as plt
import csv

class Temp_distribute1D():
    def __init__(self,dt,dx):
        self.ws = GamsWorkspace(working_directory = os.getcwd())
        self.dt = dt
        self.dx = dx
        self.drawing = 1
        
        self.jobGAMSfile = "TempDist1D_.gms"

    def solve_result(self):
        db = self.ws.add_database()
        job = self.ws.add_job_from_file(self.jobGAMSfile)
        opt = self.ws.add_options()
        opt.defines["gdxincname"] = db.name
        opt.all_model_types = "conopt"
        opt.idcgdxoutput = "test.gdx"
        dt = db.add_parameter("dt", 0)
        dt.add_record().value = self.dt
        dx = db.add_parameter("dx", 0)
        dx.add_record().value = self.dx
        job.run(opt,databases = db)
        FId = job.out_db.get_variable("obj").first_record().level

        return FId

if __name__ == '__main__' :
    test = Temp_distribute1D(0.01,0.2)
    result = test.solve_result()

    print(result)

    exit(0)