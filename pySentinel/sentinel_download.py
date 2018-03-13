#! /usr/bin/env python2
# -*- coding: iso-8859-1 -*-


from xml.dom import minidom
from datetime import date
from os.path import basename, splitext, dirname, isdir, sep
import subprocess
from re import match
import sys
import os
import multiprocessing
from functools import partial
import re
import glob
 
MAX_ATTEMPTS=5

class sentinel_configuration_download():
    
    def __init__ (self,input_file, logfile=None):
        self.MAX_DOWNLOADS=2
        self.input_file=input_file
        self.logfile=logfile
        # Variables common to any sentinel
        self.server = [
            'server',
            'user',
            'password']
        
        self.common = [
            'platformname',
            'footprint',            
            'ingestiondate',
            'beginPosition',
			'producttype',
            'endPosition']

        # Variables for Sentinel-1
        self.S1 = [
            'polarisationmode',
            'sensoroperationalmode']
        
        # Variables for Sentinel-2
        self.S2 = ['cloudcoverpercentage']

        # Variables for Sentinel-3
        self.S3 = [
           'instrumentshortname',
           'productlevel',
           'timeliness']
        
        self.additional_filters =['filename','orbitdirection']
    
    def parse_input_config(self):
        ''' Parses the information contained in a configuration file into a dictionary'''
        self.working_dir=dirname(self.input_file)
        self.config_name=splitext(basename(self.input_file))[0]
        # Prepare a list of expected input variables
        input_vars = list(self.common)
        input_vars.extend(self.server)
        input_vars.extend(self.S1)
        input_vars.extend(self.S2)
        input_vars.extend(self.S3)
        input_vars.extend(self.additional_filters)
        
        # Read contents of the configuration file
        self.config_data = dict()
        try:
            with open (self.input_file, 'r') as fid:
                for line in fid:
                    if match('\s', line):  # skip empty line
                        continue
                    elif match('#', line):  # skip comment line
                        continue
                    elif '=' in line:
                        # Remove comments in case they exist
                        line = line.split('#')[0].rstrip(' \r\n')
                        field, value = line.split('=')
                        if field in input_vars:
                            self.config_data[field] = value
        except IOError:
            print('Error reading ' + self.input_file + ' file')

        return self.config_data
    
    def get_query_command(self, query_file=None):
        '''Parses the parameters in a configuration file to a wget command line'''
        if not query_file:
            self.query_file='%s/query_results/query_%s_%s.xml'%(self.working_dir,self.config_name,date.today().strftime('%Y%m%d'))
        else:
            self.query_file=query_file
        if not isdir(self.working_dir+sep+'query_results'):
            os.mkdir(self.working_dir+sep+'query_results')
        query_command='wget --no-check-certificate --output-document=%s --user=%s --password=%s "%s/search?q='%(
                                                                 self.query_file,
                                                                 self.config_data['user'],
                                                                 self.config_data['password'],
                                                                 self.config_data['server'])
        for item in self.common:
            if item in self.config_data.keys():
                query_command+=' %s:%s AND'%(item,self.config_data[item])

        # Fill wget query depending on platform
        if self.config_data['platformname']=='Sentinel-1':
            for item in self.S1:
                if item in self.config_data.keys():
                     query_command+=' %s:%s AND'%(item,self.config_data[item])
        elif self.config_data['platformname']=='Sentinel-2':
            for item in self.S2:
                if item in self.config_data.keys():
                     query_command+=' %s:%s AND'%(item,self.config_data[item])
        elif self.config_data['platformname']=='Sentinel-3':
            for item in self.S3:
                if item in self.config_data.keys():
                     query_command+=' %s:%s AND'%(item,self.config_data[item])
        else:
            raise ValueError('wrong mission name, type Sentinel-1, Sentinel-2 or Sentinel-3')
        query_command=query_command.rstrip('AND')
        query_command+='&rows=100"'
        #print(query_command)
        proc=subprocess.Popen(query_command,shell=True,stdout=subprocess.PIPE,
                              stdin=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
        proc.wait()
    
    def run(self):
        self.parse_input_config()
        self.get_query_command()
        wget_opts={'user':self.config_data['user'],'password':self.config_data['password']}
        output_dir=self.working_dir+sep+ self.config_data['platformname']
        if self.config_data['platformname']=='Sentinel-3':
            output_dir+=sep+ self.config_data['instrumentshortname']
        
        downloaded_files=set(glob.glob(output_dir+sep+'/*.zip'))
        self.download_query_products(output_dir,
                                     wget_opts,
                                     logfile = self.logfile,
                                     downloaded_files = downloaded_files)

    def download_query_products(self, 
                                output_dir, 
                                wget_opts, 
                                logfile=None, 
                                downloaded_files=[]):
        
        jobs=[]
        download_list=[]
        
        # Get query results
        links, filenames, orbits = parse_sentinel_hub_query(self.query_file)
        
        # Loop query results and downloaded files meeting the criteria
        for link, filename, orbit in zip(links, filenames, orbits):
            output=output_dir+sep+filename
            if output+'.zip' not in downloaded_files:
                download=True
                # Filter images according to additional_filters
                if 'orbitdirection' in self.config_data.keys():
                    if orbit.lower()!=self.config_data['orbitdirection'].lower():
                        download=False
                # Filter images according to additional_filters
                if 'filename' in self.config_data.keys():
                    if not re.search(self.config_data['filename'], filename):
                        download=False
                if download ==True:
                    download_list.append(filename) 
                    jobs.append((link,output))
        
        if len(jobs)>0:
            pool = multiprocessing.Pool(processes=self.MAX_DOWNLOADS) # how much parallelism?
            pool.map(partial(wget_download_star,options=wget_opts), jobs)     
            pool.close()
            
        # Retrieve file md5 checksum from the repository    
        if logfile:
            with open(logfile,'a') as logfid:
                for download in download_list:
                    logfid.write('%s\t%s\n'%(download,True))
        
        return download_list
    
def parse_sentinel_hub_query(query_file):
    xml=minidom.parse(query_file)
    products=xml.getElementsByTagName('entry')
    
    links, filenames, orbits= [], [], []
    for prod in products:# Loop along found products
        # Get download link
        link=prod.getElementsByTagName('link')[0].attributes.items()[0][1] 
        links.append(link)
        # Find product name
        for node in prod.getElementsByTagName('str'):
            if len(node.attributes.items())>0:
                (name,field)=node.attributes.items()[0]
                if field=='filename':
                    filename= str(node.toxml()).split('>')[1].split('<')[0]
                    filenames.append(filename)
                if field=='orbitdirection':
                    orbit=str(node.toxml()).split('>')[1].split('<')[0]
                    orbits.append(orbit)
    
    return links, filenames, orbits
        
def wget_download_star(url_out,options):
    wget_download(*url_out,options_dict=options)
    return True
    
                
def wget_download(url,output,options_dict=dict()):
    if sys.platform == 'win32':
        value='$value'
    else:
        value='\$value'
    wget_command=r'wget --no-check-certificate'
    for item in options_dict.keys():
        wget_command+=' --%s=%s'%(item,options_dict[item])
    # download checksum
    url_md5=url.replace('$value','Checksum/Value/%s'%(value))
    url=url.replace('$value',value)
    wget_checksum=wget_command+' --output-document=%s.md5 "%s"'%(output,url_md5)
    proc=subprocess.Popen(wget_checksum,shell=True,stdout=subprocess.PIPE,
                               stdin=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
    proc.wait()
    with open(output+'.md5', 'r') as md5fid:
        md5=md5fid.readline()
    os.remove(output+'.md5')
    # Try to download the product to a maximum number of attemps
    attempt=1
    while attempt <= MAX_ATTEMPTS:
        wget_file=wget_command+' --output-document=%s.zip "%s"'%(output,url)
        print('Downloading file %s'%(output))
        print(wget_file)
        proc=subprocess.Popen(wget_file,shell=True,stdout=subprocess.PIPE,
                               stdin=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
        
        for line in iter(proc.stdout.readline, ''):
            print(line.rstrip('\r\n'))
        proc.stdout.close()
        proc.wait()

        # Checksum the integrity of the downloaded file
        if md5_checksum(output+'.zip', md5):
            print('Successfully downloaded %s'%(output))
            break # Download next file
        else: # Try to download again
            os.remove(output+'.zip')
            attempt+=1
    if os.path.exists(output+'.zip'):
        return True
    else:
        return False
    
def md5_checksum(infile, md5):
    
    md5_to_check=get_md5(infile).upper()
    if md5_to_check==md5:
        return True
    else:
        return False

def get_md5(fname):
    import hashlib
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()
