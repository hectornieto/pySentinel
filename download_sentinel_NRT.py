# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 09:20:01 2017

@author: hnieto
"""

from pySentinel.sentinel_download import sentinel_configuration_download
import os
import os.path as pth

if __name__=='__main__':

    workdir=os.getcwd()
    test=sentinel_configuration_download(pth.join(workdir,'Lleida_S2.txt'),logfile=pth.join(workdir,'Lleida_Download_S2A.log'))
    test.run()
    test=sentinel_configuration_download(pth.join(workdir,'Lleida_S3_SLSTR.txt'),logfile=pth.join(workdir,'Lleida_Download_S3_SLSTR.log'))
    test.run()
    test=sentinel_configuration_download(pth.join(workdir,'Lleida_S3_OLCI.txt'),logfile=pth.join(workdir,'Lleida_Download_S3_OLCI.log'))
    test.run()