'''
Created on 7 abr. 2019

@author: luis98
'''
import datetime
import os
from distutils.dir_util import copy_tree
import pytest
from connections.FTP_connection import FTP_connection


def copy_files(connection):
    direction = os.getcwd() + os.sep
    origin = direction + 'structures'
    destination = connection.anon_root + '/structures'
    copy_tree(origin, destination)


def test_0_get_all_emds_id(ftpserver):
    ftp_connection = FTP_connection("localhost", ftpserver.server_port)
    copy_files(ftpserver)
    ftp_connection.init_connection("structures")
    result = ftp_connection.get_all_emds_id("0001", "inf")
    assert result == ["0001"]


def test_1_get_all_emds_id(ftpserver):
    ftp_connection = FTP_connection("localhost", ftpserver.server_port)
    copy_files(ftpserver)
    ftp_connection.init_connection("structures")
    result = ftp_connection.get_emds_higher_than_date(
        datetime.datetime(2000, 2, 14), "0001", "inf")
    assert result == ["0001"]


def test_2_get_all_emds_id(ftpserver):
    ftp_connection = FTP_connection("localhost", ftpserver.server_port)
    copy_files(ftpserver)
    ftp_connection.init_connection("structures")
    result = ftp_connection.get_emds_higher_than_date(
        datetime.datetime(2020, 2, 14), "0001", "inf")
    assert result == []
