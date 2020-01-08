# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 08:59:33 2019

@author: TV255140
"""

from math import *
import numpy as np
from datetime import datetime, timedelta
now = datetime.now()

DeviceType      = 'TLM' # 'TLM' = Transmission Line Measurement

#%% Device dimensions
# TLM dimensions
WMesaList       = [1000,2000,5000,10000]        # All lengths in nanometers
MesaSpacing     = [0,1000000]                   # Spacing between devices when iterating through the WMesaList
LContactList    = [50,100,200,500]              # 4 contact lengths per chip
ContactSpacing  = [1000000,0]                   # Spacing between devices when iterating through the LContactList
LChannelList    = [50,100,150,200,500,1000,2000]# 8 contacts on one TLM bar, 7 steps
LEnd            = 1000                          # Length that edges of Si/Ge stick out left and right
LArm            = 750                           # Vertical length that arms stick out before bending to contacts
WArm            = 50000                         # Width of arm at the contact
Margin          = 50                            # Spill of mask 2 around the mesa

# Contact pads in a 3x3 grid
PadWidth        = 200000
PadPositions    = [[[-300000,WArm],[-300000,-WArm]],                            # Left
                   [[-300000-WArm,-300000+WArm],[-300000+WArm,-300000-WArm]],   # Bottom left
                   [[-WArm,-300000],[WArm,-300000]],                            # Bottom.. etc
                   [[300000-WArm,-300000-WArm],[300000+WArm,-300000+WArm]],     
                   [[300000,-WArm],[300000,+WArm]],
                   [[300000+WArm,300000-WArm],[300000-WArm,300000+WArm]],
                   [[WArm,300000],[-WArm,300000]],
                   [[-300000+WArm,300000+WArm],[-300000-WArm,300000-WArm]]]
ContactPadList  = [[-300000,0],[-300000,-300000],[0,-300000],[300000,-300000],[300000,0],[300000,300000],[0,300000],[-300000,300000]]

#%% Function definitions
def CurrentPosition(_StartingPoint,_StepList):
    StepArray   = np.array(_StepList)
    return [_StartingPoint[0]+sum(StepArray[:,0]),_StartingPoint[1]+sum(StepArray[:,1])]

def PosDiff(_PointA,_PointB):
    return [_PointB[0]-_PointA[0],_PointB[1]-_PointA[1]]
    
def PosSum(_PointA,_PointB):
    return [_PointB[0]+_PointA[0],_PointB[1]+_PointA[1]]

def PositionsToSteps(_PositionList):
    _StepList       = []
    for PositionIterator in range(1,len(_PositionList),1):
        _StepList   += [PosDiff(_PositionList[PositionIterator-1],_PositionList[PositionIterator])]
    return _StepList

#%% Mask 1: Si/Ge/Al, etch down Al and Ge. DEVICES
PolygonFileName = 'Polygon_%d%02d%02d_%02dh%02d_%s_Mask1.txt' % (now.year, now.month, now.day, now.hour, now.minute, DeviceType)
File            = open(PolygonFileName,"w+")
PointCounter    = 0
# Empty  list of points to contact the pads to (ends of the arms sticking out bottom and top)
ArmContactList  = [[[[[0,0],[0,0]] for ChannelIterator in range(0,len(LChannelList)+1,1)] for ContactIterator in range(0,len(LContactList),1)] for MesaIterator in range(0,len(WMesaList),1)]
for MesaIterator in range(0,len(WMesaList),1):
    for ContactIterator in range(0,len(LContactList),1):
        TotalLength     = 2*LEnd+8*LContactList[ContactIterator]+sum(LChannelList)
        DeviceCenter    = [ContactIterator*ContactSpacing[0]+MesaIterator*MesaSpacing[0],ContactIterator*ContactSpacing[1]+MesaIterator*MesaSpacing[1]]
        StartingPoint   = [ContactIterator*ContactSpacing[0]+MesaIterator*MesaSpacing[0]-TotalLength/2,ContactIterator*ContactSpacing[1]+MesaIterator*MesaSpacing[1]+WMesaList[MesaIterator]/2]
        
        # Edge that sticks out on the left
        StepList        = [[0,-WMesaList[MesaIterator]],
                            [LEnd,0],
                            [0,-LArm]]
        StartOfTrace    = CurrentPosition(StartingPoint,StepList)
        
        # Calculation of e-beam part of trace to the pads
        PadPosition1    = PosSum(DeviceCenter,PadPositions[0][0])
        PadPosition2    = PosSum(DeviceCenter,PadPositions[0][1])
        EndOfTrace      = [StartOfTrace[0]+LContactList[ContactIterator],StartOfTrace[1]-LContactList[ContactIterator]]
        Trace1          = PosSum(StartOfTrace,[x*0.04 for x in PosDiff(StartOfTrace,PadPosition1)])
        Trace2          = PosSum(EndOfTrace,[x*0.04 for x in PosDiff(EndOfTrace,PadPosition2)])
        Trace3          = PosSum(StartOfTrace,[x*0.02 for x in PosDiff(StartOfTrace,PadPosition1)])
        Trace4          = PosSum(EndOfTrace,[x*0.02 for x in PosDiff(EndOfTrace,PadPosition2)])
        ArmContactList[MesaIterator][ContactIterator][0]=[Trace3,Trace4]        
        
        StepList        += [PosDiff(StartOfTrace,Trace1),PosDiff(Trace1,Trace2),PosDiff(Trace2,EndOfTrace)]
        StepList        += [[0,LArm+LContactList[ContactIterator]]]
        
        # Arms below
        for ChannelIterator in range(0,len(LChannelList),2):
            if len(LChannelList)%2 or ChannelIterator<len(LChannelList)-1:
                StepList    += [[LChannelList[ChannelIterator],0],
                                 [LContactList[ContactIterator],0]]
            if ChannelIterator<len(LChannelList)-2:
                StepList.append([LChannelList[ChannelIterator+1],0])
                StepList    += [[0,-(1+(ChannelIterator==0 or ChannelIterator==2))*LArm]] # Long arms in the middle
                StartOfTrace= CurrentPosition(StartingPoint,StepList)
                
                # Calculation of e-beam part of trace to the pads
                PadPosition1= PosSum(DeviceCenter,PadPositions[1+ChannelIterator/2][0])
                PadPosition2= PosSum(DeviceCenter,PadPositions[1+ChannelIterator/2][1])
                EndOfTrace  = [StartOfTrace[0]+LContactList[ContactIterator],StartOfTrace[1]-(1-ChannelIterator/2)*LContactList[ContactIterator]]
                Trace1      = PosSum(StartOfTrace,[x*0.04 for x in PosDiff(StartOfTrace,PadPosition1)])
                Trace2      = PosSum(EndOfTrace,[x*0.04 for x in PosDiff(EndOfTrace,PadPosition2)])
                Trace3      = PosSum(StartOfTrace,[x*0.02 for x in PosDiff(StartOfTrace,PadPosition1)])
                Trace4      = PosSum(EndOfTrace,[x*0.02 for x in PosDiff(EndOfTrace,PadPosition2)])
                ArmContactList[MesaIterator][ContactIterator][1+ChannelIterator/2]=[Trace3,Trace4]
                
                StepList    += [PosDiff(StartOfTrace,Trace1),PosDiff(Trace1,Trace2),PosDiff(Trace2,EndOfTrace)]
                StepList    += [[0,(1+(ChannelIterator==0 or ChannelIterator==2))*LArm+(1-ChannelIterator/2)*LContactList[ContactIterator]]]  # Long arms in the middle
        
        # Edge that sticks out on the right
        StepList            += [[LEnd,0],
                                [0,WMesaList[MesaIterator]],
                                [-LEnd,0]]
        if (len(LChannelList)+1)%2:     # If there is an odd number of arms, take a step first
            StepList.append([-LChannelList[len(LChannelList)-1]-LContactList[ContactIterator],0])
            
        # Arms above
        for ChannelIterator in range(len(LChannelList)-len(LChannelList)%2,-1,-2):
            StepList        += [[0,(1+(ChannelIterator==2 or ChannelIterator==4))*LArm]]    # Long arms in the middle
            StartOfTrace    = CurrentPosition(StartingPoint,StepList)
            
            # Calculation of e-beam part of trace to the pads
            PadPosition1    = PosSum(DeviceCenter,PadPositions[7-ChannelIterator/2][0])
            PadPosition2    = PosSum(DeviceCenter,PadPositions[7-ChannelIterator/2][1])
            EndOfTrace      = [StartOfTrace[0]-LContactList[ContactIterator],StartOfTrace[1]+(ChannelIterator/2-1-(ChannelIterator==6))*LContactList[ContactIterator]]
            Trace1          = PosSum(StartOfTrace,[x*0.04 for x in PosDiff(StartOfTrace,PadPosition1)])
            Trace2          = PosSum(EndOfTrace,[x*0.04 for x in PosDiff(EndOfTrace,PadPosition2)])
            Trace3          = PosSum(StartOfTrace,[x*0.02 for x in PosDiff(StartOfTrace,PadPosition1)])
            Trace4          = PosSum(EndOfTrace,[x*0.02 for x in PosDiff(EndOfTrace,PadPosition2)])
            ArmContactList[MesaIterator][ContactIterator][7-ChannelIterator/2]=[Trace3,Trace4]        
            
            StepList        += [PosDiff(StartOfTrace,Trace1),PosDiff(Trace1,Trace2),PosDiff(Trace2,EndOfTrace)]
            StepList        += [[0,-(1+(ChannelIterator==2 or ChannelIterator==4))*LArm-(ChannelIterator/2-1-(ChannelIterator==6))*LContactList[ContactIterator]],
                                 [-LChannelList[ChannelIterator],0],
                                 [-LContactList[ContactIterator],0]]
            if ChannelIterator>1:
                StepList.append([-LChannelList[ChannelIterator-1],0])
        StepList.append([-LEnd,0])
                                
        # Convert list of steps to list of points and write them to the file
        PointList       = [StartingPoint]
        for Step in StepList:
            PointList.append([PointList[len(PointList)-1][0]+Step[0],PointList[len(PointList)-1][1]+Step[1]])
        for Point in PointList:
            File.write("%d	%d\n" % (Point[0],Point[1]))
            PointCounter+=1
        File.write("\n")
    
    # Go back to the start of the row
    TotalLength         = 2*LEnd+8*LContactList[0]+sum(LChannelList)
    StartingPoint       = [-TotalLength/2,ContactIterator*ContactSpacing[1]+MesaIterator*MesaSpacing[1]+WMesaList[MesaIterator]/2]
    File.write("%d	%d\n" % (StartingPoint[0],StartingPoint[1]))
    PointCounter+=1
    File.write("\n")
File.close()
print('File:\n    '+PolygonFileName+'\n')
print("Total number of points:\n    %d" % PointCounter)

#%% Mask 1: Si/Ge/Al, etch down Al and Ge. CONTACT PADS
PolygonFileName = 'Polygon_%d%02d%02d_%02dh%02d_%s_Mask1_CONTACT_PADS.txt' % (now.year, now.month, now.day, now.hour, now.minute, DeviceType)
File            = open(PolygonFileName,"w+")
PointCounter    = 0
# Empty  list of points to contact the pads to (ends of the arms sticking out bottom and top)
for MesaIterator in range(0,len(WMesaList),1):
    for ContactIterator in range(0,len(LContactList),1):
        TotalLength     = 2*LEnd+8*LContactList[ContactIterator]+sum(LChannelList)
        DeviceCenter    = [ContactIterator*ContactSpacing[0]+MesaIterator*MesaSpacing[0],ContactIterator*ContactSpacing[1]+MesaIterator*MesaSpacing[1]]
        StartingPoint   = [ContactIterator*ContactSpacing[0]+MesaIterator*MesaSpacing[0]-TotalLength/2,ContactIterator*ContactSpacing[1]+MesaIterator*MesaSpacing[1]+WMesaList[MesaIterator]/2]
        StepList        = []
        
        # Contact pads
        for PadCenter in ContactPadList:
            StepList    += PositionsToSteps([StartingPoint,
                                PosSum(PosSum(DeviceCenter,PadCenter),[PadWidth/2,PadWidth/2]),
                                PosSum(PosSum(DeviceCenter,PadCenter),[-PadWidth/2,PadWidth/2]),
                                PosSum(PosSum(DeviceCenter,PadCenter),[-PadWidth/2,-PadWidth/2]),
                                PosSum(PosSum(DeviceCenter,PadCenter),[PadWidth/2,-PadWidth/2]),
                                PosSum(PosSum(DeviceCenter,PadCenter),[PadWidth/2,PadWidth/2]),
                                StartingPoint])
        for PadIterator in range(0,len(LChannelList)+1,1):
            StepList    += PositionsToSteps([StartingPoint,
                                ArmContactList[MesaIterator][ContactIterator][PadIterator][0],
                                PosSum(DeviceCenter,PadPositions[PadIterator][0]),
                                PosSum(DeviceCenter,PadPositions[PadIterator][1]),
                                ArmContactList[MesaIterator][ContactIterator][PadIterator][1],
                                ArmContactList[MesaIterator][ContactIterator][PadIterator][0],
                                StartingPoint])
                                
        # Convert list of steps to list of points and write them to the file
        PointList       = [StartingPoint]
        for Step in StepList:
            PointList.append([PointList[len(PointList)-1][0]+Step[0],PointList[len(PointList)-1][1]+Step[1]])
        for Point in PointList:
            File.write("%d	%d\n" % (Point[0],Point[1]))
            PointCounter+=1
        File.write("\n")
    
    # Go back to the start of the row
    TotalLength         = 2*LEnd+8*LContactList[0]+sum(LChannelList)
    StartingPoint       = [-TotalLength/2,ContactIterator*ContactSpacing[1]+MesaIterator*MesaSpacing[1]+WMesaList[MesaIterator]/2]
    File.write("%d	%d\n" % (StartingPoint[0],StartingPoint[1]))
    PointCounter+=1
    File.write("\n")
File.close()
print('File:\n    '+PolygonFileName+'\n')
print("Total number of points:\n    %d" % PointCounter)

#%% Mask 2: Si/Ge/Al, etch down only Al.
PolygonFileName = 'Polygon_%d%02d%02d_%02dh%02d_%s_Mask2.txt' % (now.year, now.month, now.day, now.hour, now.minute, DeviceType)
File            = open(PolygonFileName,"w+")
PointCounter    = 0
for MesaIterator in range(0,len(WMesaList),1):
    for ContactIterator in range(0,len(LContactList),1):
        TotalLength     = 2*LEnd+2*Margin+8*LContactList[ContactIterator]+sum(LChannelList)
        DeviceCenter    = [ContactIterator*ContactSpacing[0]+MesaIterator*MesaSpacing[0],ContactIterator*ContactSpacing[1]+MesaIterator*MesaSpacing[1]]
        StartingPoint   = [ContactIterator*ContactSpacing[0]+MesaIterator*MesaSpacing[0]-TotalLength/2,ContactIterator*ContactSpacing[1]+MesaIterator*MesaSpacing[1]+WMesaList[MesaIterator]/2+Margin]
        
        # Edge that sticks out on the left
        StepList        = [[0,-WMesaList[MesaIterator]-2*Margin],
                            [LEnd+Margin,0],
                            [0,WMesaList[MesaIterator]+2*Margin],
                            [-LEnd-Margin,0],
                            [LEnd+Margin+LContactList[ContactIterator],0]]
        
        # Arms below
        for ChannelIterator in range(0,len(LChannelList),1):
            StepList    += [[0,-WMesaList[MesaIterator]-2*Margin],
                            [LChannelList[ChannelIterator],0],
                            [0,WMesaList[MesaIterator]+2*Margin],
                            [-LChannelList[ChannelIterator],0],
                            [LChannelList[ChannelIterator]+LContactList[ContactIterator],0]]
        
        # Edge that sticks out on the right
        StepList        += [[0,-WMesaList[MesaIterator]-2*Margin],
                            [LEnd+Margin,0],
                            [0,WMesaList[MesaIterator]+2*Margin],
                            [-LEnd-Margin,0],
                            [LEnd+Margin-TotalLength,0]]
            
        # Convert list of steps to list of points and write them to the file
        PointList       = [StartingPoint]
        for Step in StepList:
            PointList.append([PointList[len(PointList)-1][0]+Step[0],PointList[len(PointList)-1][1]+Step[1]])
        for Point in PointList:
            File.write("%d	%d\n" % (Point[0],Point[1]))
            PointCounter+=1
        File.write("\n")
    
    # Go back to the start of the row
    TotalLength     = 2*LEnd+2*Margin+8*LContactList[0]+sum(LChannelList)
    StartingPoint   = [-TotalLength/2,ContactIterator*ContactSpacing[1]+MesaIterator*MesaSpacing[1]+WMesaList[MesaIterator]/2+Margin]
    File.write("%d	%d\n" % (StartingPoint[0],StartingPoint[1]))
    PointCounter+=1
    File.write("\n")
File.close()
print('File:\n    '+PolygonFileName+'\n')
print("Total number of points:\n    %d" % PointCounter)
