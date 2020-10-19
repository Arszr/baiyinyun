#! /usr/bin/python
#coding: utf-8
import os
import sys
import zipfile
 
path = 'e:\\pythonWorkSpace\\practice'
zipfilename = 'e:\\version.zip'
 
# path is a direactory or not.
if not os.path.isdir(path):
    print path + ' No such a direactory'
    exit()
 
if os.path.exists(zipfilename):
    # zipfilename is exist.Append.
    print 'Add files into ' + zipfilename
    zipfp = zipfile.ZipFile(zipfilename, 'a' ,zipfile.ZIP_DEFLATED)
    for dirpath, dirnames, filenames in os.walk(path, True):
        baifor filaname in filenames:
            direactory = os.path.join(dirpath,filaname)
            print 'Add... ' + direactory
            zipfp.write(direactory)
else:
    # zipfilename is not exist.Create.
    print 'Create new file ' + zipfilename
    zipfp = zipfile.ZipFile(zipfilename, 'w' ,zipfile.ZIP_DEFLATED)
    for dirpath, dirnames, filenames in os.walk(path, True):
        for filaname in filenames:
            direactory = os.path.join(dirpath,filaname)
            print 'Compress... ' + direactory
            zipfp.write(direactory)
# Flush and Close zipfilename at last. 
zipfp.close()