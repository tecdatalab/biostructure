'''
Created on 30 mar. 2019

@author: luis98
'''
from ftplib import FTP
import datetime

def max_datetime(date1, date2):
    if date1 == None:
        return date2
    elif date2 == None:
        return date1
    elif date1 > date2:
        return date1
    else:
        return date2

class FTP_connection(object):
    '''
    classdocs
    '''
    __ftp = None

    def __init__(self,ftp_server = 'ftp.ebi.ac.uk', port = None):
        if port == None:
            self.__ftp = FTP()  
            self.__ftp.connect(ftp_server)
        else:
            self.__ftp = FTP()  
            self.__ftp.connect(ftp_server, port)
    
    def init_connection(self,initial_directory='/pub/databases/emdb/structures/'):
        self.__ftp.login()
        if initial_directory != None:
            self.__ftp.cwd(initial_directory)  # change into "structures" directory     
    
    
    def close_connection(self):
        self.__ftp.close()
        
        
    def get_all_emds_id(self):
        files = []
        result = []
        self.__ftp.retrlines("NLST",files.append)
        for i in range(len(files)):
            try:
                self.__ftp.voidcmd("MDTM EMD-{0}/map/emd_{0}.map.gz".format(files[i].replace("EMD-", "")))
                result.append(files[i].replace("EMD-", ""))
            except:
                pass
            
            return result
    
    
    def get_emds_higher_than_date(self,date):
        emds_id = self.get_all_emds_id()
        result = []
        for emd in emds_id:
            time_xml = None
            time_image = None
            time_map = None
            try:
                time_xml = self.__ftp.voidcmd("MDTM EMD-{0}/header/emd-{0}.xml".format(emd))
                time_xml = datetime.datetime.strptime(time_xml[4:], "%Y%m%d%H%M%S")
            except:
                pass
            try:
                time_image = self.__ftp.voidcmd("MDTM EMD-{0}/images/emd_{0}.png".format(emd))
                time_image = datetime.datetime.strptime(time_image[4:], "%Y%m%d%H%M%S")
            except:
                pass
            try:
                time_map = self.__ftp.voidcmd("MDTM EMD-{0}/map/emd_{0}.map.gz".format(emd))
                time_map = datetime.datetime.strptime(time_map[4:], "%Y%m%d%H%M%S")
            except:
                pass
                
            max_value = max_datetime(max_datetime(time_xml,time_image),time_map)
            
            if max_value > date:
                result.append(emd)
        
        return result
    

    def get_ftp(self):
        return self.__ftp


    def set_ftp(self, value):
        self.__ftp = value


    def del_ftp(self):
        del self.__ftp


    ftp = property(get_ftp, set_ftp, del_ftp, "ftp's docstring")

'''
test = FTP_connection()
test.init_connection()
#directorios = test.get_emds_higher_than_date(datetime.datetime(2019, 6, 1))
directorios = test.get_all_emds_id()
print( len(directorios) )
print(directorios[0])
test.close_connection()
'''
