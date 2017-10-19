import threading
import matlab.engine
import time
import datetime


class matlab_start(threading.Thread):

    eg = matlab.engine

    def __init__(self):
        threading.Thread.__init__(self)

    def run(self):
        matlab_start.eg = matlab_start.eg.start_matlab()
