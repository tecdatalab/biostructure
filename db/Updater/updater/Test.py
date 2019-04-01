'''
Created on 31 mar. 2019

@author: luis98
'''
import requests

URL = "http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0001/header/em4d-0002.xml"

response = requests.get(URL)
with open('feed.xml', 'wb') as file:
    file.write(response.content)