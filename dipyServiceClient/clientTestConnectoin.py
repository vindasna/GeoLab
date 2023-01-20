#!/usr/bin/python3
import sys, os, subprocess, shutil

import time


import socket



def client_program() :
    ################################## Client ##################################
    host = socket.gethostname()  # as both code is running on same pc
    port = 5000  # socket server port number

    client_socket = socket.socket()  # instantiate
    client_socket.connect( ( host, port ) )  # connect to the server

    maxBytesReceive = 1024
    closeClient = False
    client_socket.send( "test".encode() )  # send message
    while not closeClient :
        data = client_socket.recv( maxBytesReceive ).decode()  # receive response
        data = data.split( " : " )
        if ( data[ 0 ] == "0" ) :
            closeClient = True
            print( "Connection to server-client : OK" )
            client_socket.close()  # close the connection
        elif ( data[ 0 ] == "Number client connected" ) :
            print( f"Number of clients connected to server {data[ 1 ]}" )
        elif ( data[ 0 ] == "Number total connections" ) :
            print( f"Number of total connections to server {data[ 1 ]}" )


if __name__ == '__main__':
    t1 = time.time()
    client_program()
    t2 = time.time()
    print( f"Time ellapsed : {t2 - t1}" )
