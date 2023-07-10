#!/usr/bin/env python
# -*- encoding: utf-8 -*-

"""
@File    :   EncoderMikeController.py
@Time    :   2022/08/09 11:59:58
@Author  :   Roney D. Silva
@Contact :   roneyddasilva@gmail.com
"""
from __future__ import absolute_import
import serial
import re
import time
import numpy as np
from collections import namedtuple
import serial.tools.list_ports

# Global Variables
## Namedtuple for setup of serial conections. PID values were found with tests and the speed speed of the settings on each device
Port = namedtuple("Port", ["port", "pid", "speed"])
## Tuple with Encoder() serial info.
PORT_ENCODER = Port("", 8963, 9600)
## Tuple with Dynamometer() serial info.
PORT_DYNAMOMETER = Port("", 24577, 38400)


## The findDevices function is used to find the USB ports that are connected to the encoder and dynamometer.
# It checks if they are connected by checking their PID, which is unique for each device.
# If it finds a port with a matching PID, it will update its corresponding global variable with the new port name.
# @return
#     A tuple with the ports that are connected to a device
def findDevices():
    global PORT_ENCODER
    global PORT_DYNAMOMETER
    listOfConnected = list(serial.tools.list_ports.comports())
    for p in listOfConnected:
        # check if Enconder it is connected and which is this USB number
        if p.pid == PORT_ENCODER.pid:
            PORT_ENCODER = PORT_ENCODER._replace(port=p.device)
            print("Encoder Found in", p.device)
        elif p.pid == PORT_DYNAMOMETER.pid:
            PORT_DYNAMOMETER = PORT_DYNAMOMETER._replace(port=p.device)
            print("Dynamometer Found in", p.device)
        else:
            print("Do not found any device on USB")


findDevices()


class EncoderMikeController:

    # Class variables
    ## Variable to get all messagens on serial channel.
    bufferMessage = ""
    ## Index of Home
    indexOfHome = 0.0
    ## Index for set Zero Strain
    indexOfZeroStrain = 0
    ## List with absolute position
    absolutePosition = [0,0,0]
    ## List with strain
    strain = [0,0,0]
    ## Initial traction element length.
    initialLength = 0.0
    ## Section area
    sectionArea = 0.0
    ## \todo: Implement a methods to find indexOfZeroStrain in a safe way.
    
    ## __init__ The constructor
    # @@param speed (int, optional): Baudrate to use. Defaults to 9600.
    #     sizeInterval(int, optional): Distânce between indexOfZeroStrain in micrometers.
    #     indexOfZeroStrain(float, optional): Initial value to indexOfZeroStrain.
    #
    def __init__(self, initialLength, indexOfZeroStrain, speed=9600):
        self.initialLength = float(initialLength)
        self.indexOfZeroStrain=indexOfZeroStrain
        self.indexOfHome = indexOfZeroStrain+initialLength
        # Change the speed os serial if necessary.
        self.PORT_ENCODER = PORT_ENCODER._replace(speed=speed)
        self.ser = serial.Serial(
            self.PORT_ENCODER.port,
            self.PORT_ENCODER.speed,
            parity=serial.PARITY_NONE,
            stopbits=serial.STOPBITS_ONE,
            bytesize=serial.EIGHTBITS,
        )
        # time.sleep(2)
        self.setPositionMode()
        # time.sleep(2)

    ## Set position mode for all axes
    def setPositionMode(self):
        self.sendCommand(1, "PM")
        self.sendCommand(2, "PM")
        self.sendCommand(3, "PM")

    ## """
    # The setVelocityMode function sets the velocity on all axes.
    # @param
    #     self: Refer to the object instance itself, and is used to access variables that belongs to the class
    def setVelocityMode(self):
        # set Velocity mode for each axis.
        self.sendCommand(1, "VM")
        self.sendCommand(2, "VM")
        self.sendCommand(3, "VM")

    ## sendCommand Function for sending command to controller.
    # @param
    #     axis (Interger,float): 1 for X, 2 for Y and 3 for Z.
    #     word (String): Commands that must be consisted in the manual
    #     n (str, optional): The numeric value to command. For example, a certain desired speed or position value . Defaults to "".
    def sendCommand(self, axis, word, n=""):
        self.ser.write((str(axis) + word.upper() + str(n) + "\r").encode())
        # It awaits 100 milliseconds
        time.sleep(0.1)

    # """
    # getInformation Function to collect information such as position and speed
    # @param
    #     axis (Interger): 1 for X, 2 for Y and 3 for Z.
    #     word (String): Commands that must be consisted in the manual.
    # @return
    #     Float : Tha value of requisition like position and velocity.
    # """
    def getInformation(self, axis, word):
        for tryCount in range(10):
            # send a command
            self.sendCommand(axis, word)
            # read the buffer and split the messages.
            # Use Regex to find float
            self.bufferMessage = str(self.ser.read_all())
            # Get the last value on message and convert to float
            try:
                # position = float(self.bufferMessage[-2].split(" ")[-1])
                position = float(re.findall("[-]?\d+\.\d+",str(self.bufferMessage))[-1])
                return position
            except:
                # waint for another try
                time.sleep(0.1)


    ## goToX Move the axis X to absolute position
    # @param
    #     position (interger): absolute position
    def goToX(self, position):
        self.setPositionMode()
        time.sleep(0.1)
        self.sendCommand(1, "ma", position)

    ## goToY Move  the axis Y to absolute position
    # @param
    #    position (interger): absolute position
    def goToY(self, position):
        self.setPositionMode()
        time.sleep(0.1)
        self.sendCommand(2, "ma", position)


    ## goToZ Move  the axis Z to absolute position
    # @param
    #     position (interger): absolute position
    def goToZ(self, position):
        self.setPositionMode()
        time.sleep(0.1)
        self.sendCommand(3, "ma", position)


    ## walkX Move the axis X micrometers
    # @param
    #     distance (interger): value postive or negative micrometeres
    def walkX(self, distance):
        self.setPositionMode()
        time.sleep(0.1)
        self.sendCommand(1, "mr", distance)

    ## walkY Move the axis Y in micrometers
    # @param
    #     distance (interger): value postive or negative micrometeres
    def walkY(self, distance):
        self.setPositionMode()
        time.sleep(0.1)
        self.sendCommand(2, "mr", distance)

    ## walkZ Move the axis Z in micrometers
    # @param
    #     distance (interger): value postive or negative micrometeres
    def walkZ(self, distance):
        self.setPositionMode()
        time.sleep(0.1)
        self.sendCommand(3, "mr", distance)

    
    ## The getPositionX function returns the position of the first motor.
    # @param
    #     self: Refer to the object itself
    # @return
    #     The position of the first element in the list
    def getPositionX(self):

        self.absolutePosition[0] = self.getInformation(1, "TP")
        return self.absolutePosition[0]



    ## The getPositionY function returns the position of the axis 2(Y).
    # @param
    #     self: Refer to the object itself
    # @return:
    #     The position of the  of the axis 2(Y)
    def getPositionY(self):
        self.absolutePosition[1]  = self.getInformation(2, "TP")
        return self.absolutePosition[1]

    ## The getPositionZ function returns the position of the axis 3(Z).
    # @param
    #     self: Refer to the object itself
    # @return:
    #     The position of the  of the axis 3(Z)
    def getPositionZ(self):
        self.absolutePosition[2] = self.getInformation(3, "TP")
        return self.absolutePosition[2]

    
    ## The function get_index_zero_stran() update #indexOfZeroStrain variable.
    # @param
    #       axis (Interger): 1 for X, 2 for Y and 3 for Z.
    
    def get_index_zero_stran(self,axis):
        if axis == 1:
            self.indexOfZeroStrain =  self.getPositionX()
        elif axis == 2:
            self.indexOfZeroStrain =  self.getPositionY()
        elif axis == 3:
            self.indexOfZeroStrain =  self.getPositionZ()
        
    
    ## The function goHome2() move the axis 2 to Home position.
    def goHomeY(self):
        self.goToY(self.indexOfHome)


    ## The putEncoderRelativeToZeroStrainPoint function move to #distance with respect to zero strain point
    # @param
    #     distance (float, optional): The distance to move. Defaults to 3000.
    def putEncoderRelativeToZeroStrainPoint(self, distance=3000):
        self.goToY(self.indexOfZeroStrain + distance)

    def calcStrainY(self):
        self.getPositionY()
        self.strain[1] = (self.absolutePosition[1] - self.indexOfHome)/self.initialLength
        return self.strain[1]

## \note: To avoid having to empty the buffer all the time, it was chosen to open the serial connection before readings and close it when using them.
    
## \note: To avoid the encoder mike controller change to velocity mode, we force to setPosition in each goTO ans walk command.

## \todo: Treat issues of invalid measurements more robustly.

## \todo: Implement the option to receive measurement only when requested.

## The Dynamomenter Class contains methods to configure, collect and control the dynamometer.    
class Dynamometer:
    global PORT_DYNAMOMETER
    strain_value = 0.0
    ## The constructor.
    # @param
    #    self: Refer to the object itself
    #    speed=38400: Change the default speed of the serial port.
    def __init__(self, speed=38400):
        
        ## change speed if necessary
        self.PORT_DYNAMOMETER = PORT_DYNAMOMETER._replace(speed=speed)
        self.ser = serial.Serial(
            self.PORT_DYNAMOMETER.port, self.PORT_DYNAMOMETER.speed
        )
        self.ser.close()
        time.sleep(1)
        # With the unit of measure in Newtons, the format of String Str (Ser.Readline ()). Split ("") is ["B '",' ',', 'F.FFN', '', '' ', "\\ r \\ n '"].So, to collect the measurement use Float (Reading[3][:-1])
        self.scaleAdjustForNewton()

    ## The update_strain function update the strain_value variable read from the dynamometer.
    # It is a float, but it may be None if no data was received.
    # @param self: Access variables that belongs to the class.
    #
    def update_strain(self):
        while True:
            try:
                self.ser.open()
                time.sleep(0.05)
                # self.strain_value = float(str(self.ser.read_all()).split(" ")[-5][:-1])
                # Use regex to find float number with the signal.
                self.strain_value = float(re.findall("[-+]?\d+\.\d+", str(self.ser.read_all()))[-1])
                self.ser.close()
                break
            except:# ValueError or IndexError:
                self.ser.close()

    ## The set_zero function sets  of the dynamometer.
    # @param
    #     self: Access the class attributes
    def set_zero(self):
        self.ser.open()
        time.sleep(1)
        self.ser.write("Z\r\n".encode())
        self.ser.close()

    ## The scaleAdjustForNewton function is used to adjust the scale for Newton.
    # It sends a command to the scale, waits 1 second, and then reads in all of the data from that serial port.
    # The function checks if there are 6 values in this data string (the first value is always 'S').  If there are not 6 values, it will send another command (&quot;e&quot;) and wait 1 more second before reading again.  This process repeats until either 6 values are read.
    # @param
    #     self: Access the variables and methods of the class in python
    #
    def scaleAdjustForNewton(self):
        print("Adjust scale for Newton")
        # self.ser.reset_input_buffer()
        while True:
            self.ser.open()
            self.ser.write("E\r\n".encode())
            time.sleep(.1)
            self.ser.reset_input_buffer()
            _read = len(str(self.ser.readline()).split(" "))
            if _read >= 6:
                self.ser.close()
                break
            else:
                print(_read)
                self.ser.close()


