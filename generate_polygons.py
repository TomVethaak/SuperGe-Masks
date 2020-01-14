# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 08:59:33 2019

@author: TV255140
"""

from math import *
import numpy as np
from datetime import datetime, timedelta
now = datetime.now()

DeviceType      = 'TLM_hall_bar' # 'TLM' = Transmission Line Measurement

#%% Device dimensions
# Chip dimensions
ChipWidth               = 6000000
ChipLength              = 8000000
HorizontalCells         = 6
HorizontalCellSpacing   = 1000000
VerticalCells           = 8
VerticalCellSpacing     = -1000000

# Alignment marks
GlobalCrossThickness    = 10000
GlobalCrossWidth        = 500000
GlobalCrossLongLength   = 2000000
GlobalCrossPoints       = [[-ChipWidth/2,ChipLength/2],                             [ChipWidth/2,ChipLength/2],
                           [-ChipWidth/2,0],                                        [ChipWidth/2,0],
                           [-ChipWidth/2,-ChipLength/2],    [0,-ChipLength/2],      [ChipWidth/2,-ChipLength/2],]
GlobalLongCrossPoints   = [[0,ChipLength/2]]
ChipCrossThickness      = 1000
ChipCrossWidth          = 15000
ChipMarkGroupHSpacing   = 1000000
ChipMarkGroupVSpacing   = 2000000
ChipMarkGroupPoints     = [[-ChipMarkGroupHSpacing,ChipMarkGroupVSpacing],    [ChipMarkGroupHSpacing,ChipMarkGroupVSpacing],
                                                                        [0,0],
                           [-ChipMarkGroupHSpacing,-ChipMarkGroupVSpacing],   [ChipMarkGroupHSpacing,-ChipMarkGroupVSpacing]]
ChipMarkSetSpacing      = 100000
ChipMarkSetPoints       = [[-ChipMarkSetSpacing,ChipMarkSetSpacing],    [ChipMarkSetSpacing,ChipMarkSetSpacing],
                                                                    [0,0],
                           [-ChipMarkSetSpacing,-ChipMarkSetSpacing],   [ChipMarkSetSpacing,-ChipMarkSetSpacing]]
ChipMarkSpacing         = 40000
ChipMarkPoints          = [[-ChipMarkSpacing/2,ChipMarkSpacing/2],  [ChipMarkSpacing/2,ChipMarkSpacing/2],
                           [-ChipMarkSpacing/2,-ChipMarkSpacing/2], [ChipMarkSpacing/2,-ChipMarkSpacing/2]]

# Contact pads in a 3x3 grid
PadSpacing              = 300000
PadWidth                = 200000
Litho1Fraction          = 0.04                          # Hi-res litho step: until what fraction of distance to pads
Litho2Fraction          = 0.02                          # Lo-res litho step: from what fraction of distance to pads
PadPositions            = [[[-PadSpacing,WArm],[-PadSpacing,-WArm]],                            # Left
                           [[-PadSpacing-WArm,-PadSpacing+WArm],[-PadSpacing+WArm,-PadSpacing-WArm]],   # Bottom left
                           [[-WArm,-PadSpacing],[WArm,-PadSpacing]],                            # Bottom.. etc
                           [[PadSpacing-WArm,-PadSpacing-WArm],[PadSpacing+WArm,-PadSpacing+WArm]],
                           [[PadSpacing,-WArm],[PadSpacing,+WArm]],
                           [[PadSpacing+WArm,PadSpacing-WArm],[PadSpacing-WArm,PadSpacing+WArm]],
                           [[WArm,PadSpacing],[-WArm,PadSpacing]],
                           [[-PadSpacing+WArm,PadSpacing+WArm],[-PadSpacing-WArm,PadSpacing-WArm]]]
ContactPadList          = [[-PadSpacing,0],[-PadSpacing,-PadSpacing],[0,-PadSpacing],[PadSpacing,-PadSpacing],[PadSpacing,0],[PadSpacing,PadSpacing],[0,PadSpacing],[-PadSpacing,PadSpacing]]

# TLM dimensions
TLMTopLeftCoordinate    = [0,0]
WMesaList               = [1000,2000,5000,10000]        # All lengths in nanometers
MesaSpacing             = [0,VerticalCellSpacing]       # Spacing between devices when iterating through the WMesaList
LContactList            = [50,100,200,500]              # 4 contact lengths per chip
ContactSpacing          = [HorizontalCellSpacing,0]     # Spacing between devices when iterating through the LContactList
LChannelList            = [50,100,150,200,500,1000,2000]# 8 contacts on one TLM bar, 7 steps
LEnd                    = 1000                          # Length that edges of Si/Ge stick out left and right
LArm                    = 750                           # Vertical length that arms stick out before bending to contacts
WArm                    = 50000                         # Width of arm at the contact
TLMMargin               = 100                           # Spill of mask 2 around the mesa

# Hall bars
HallTopLeftCoordinate   = [4,0]
HallMesaSpacing         = MesaSpacing
HallWMesaList           = [5000,10000,20000,50000]
HallLayerSpacing        = ContactSpacing
HallLayerList           = [1,0]
HallArmSpacing          = 50000
HallLMesa               = 200000
HallWContact            = 2000
HallLArm                = 3000
HallMargin              = 0

# Greek crosses
GreekTopLeftCoordinate  = [0,4]
GreekWidthSpacing       = MesaSpacing
GreekWidthList          = [2000,5000,10000,20000]
GreekLengthSpacing      = ContactSpacing
GreekLengthList         = [20000,20000,20000,20000,20000,50000]
GreekLayerList          = [1,0]
GreekMargin             = 0


#%% Function definitions
# Parameter calculations
def CoordinateToPosition(Coordinate):
    return [((0.5-HorizontalCells/2.)+Coordinate[0])*HorizontalCellSpacing,((Coordinate[1]-VerticalCells/2.+0.5)*VerticalCellSpacing)]

# Single point operations
def IntermediatePoint(_PointA,_PointB,_Fraction):
    return PosSum(_PointA,PosDiffFactor(_PointA,_PointB,_Fraction))

def MoveTo(_ReferencePoint,_StepList,_Point):
    return [PosDiff(StepsToPosition(_ReferencePoint,_StepList),_Point)]

def PosDiff(_PointA,_PointB):
    return [_PointB[0]-_PointA[0],_PointB[1]-_PointA[1]]

def PosDiffFactor(_PointA,_PointB,_Fraction):
    return [x*_Fraction for x in PosDiff(_PointA,_PointB)]

def PosSum(_PointA,_PointB):
    return [_PointB[0]+_PointA[0],_PointB[1]+_PointA[1]]

def PosListSum(_PosList):
    _PosArray       = np.array(_PosList)
    return [sum(_PosArray[:,0]),sum(_PosArray[:,1])]

def StepFactor(_Step,_Factor):
    return [x*_Factor for x in _Step]

def StepListFactor(_StepList,_Factor):
    return [[x*_Factor for x in _Step] for _Step in _StepList]

# List operations
def PositionsToSteps(_PositionList):
    _StepList       = []
    for PositionIterator in range(1,len(_PositionList),1):
        _StepList   += [PosDiff(_PositionList[PositionIterator-1],_PositionList[PositionIterator])]
    return _StepList
    
def StepsToPosition(_StartingPoint,_StepList):
    StepArray       = np.array(_StepList)
    return [_StartingPoint[0]+sum(StepArray[:,0]),_StartingPoint[1]+sum(StepArray[:,1])]

# Path macros
def CrossSteps(_CrossThickness,_CrossWidth):
    _CS            = _CrossWidth-_CrossThickness            # Cross Step size
    _StepList       = [[_CrossThickness/2,_CrossThickness/2]]
    for _i in range(0,4,1):
        _StepList   += [[(_i%2)*((-1)**(1+_i//2))*_CS,(1-_i%2)*((-1)**(_i//2))*_CS],
                        [(1-_i%2)*((-1)**(1+_i//2))*_CrossThickness,(_i%2)*((-1)**(1+_i//2))*_CrossThickness],
                        [(_i%2)*((-1)**(_i//2))*_CS,(1-_i%2)*((-1)**(1+_i//2))*_CS]]
    _StepList       += [[-_CrossThickness/2,-_CrossThickness/2]]
    return _StepList

def LongCrossSteps(_CrossThickness,_CrossWidth,_CrossLength):
    _GCHS           = _CrossLength-_CrossThickness            # Global Cross Horizontal Step size
    _GCVS           = _CrossWidth-_CrossThickness
    _StepList       = [[_CrossThickness/2,_CrossThickness/2]]
    for _i in range(0,4,1):
        _StepList   += [[(_i%2)*((-1)**(1+_i//2))*_GCHS,(1-_i%2)*((-1)**(_i//2))*_GCVS],
                        [(1-_i%2)*((-1)**(1+_i//2))*_CrossThickness,(_i%2)*((-1)**(1+_i//2))*_CrossThickness],
                        [(_i%2)*((-1)**(_i//2))*_GCHS,(1-_i%2)*((-1)**(1+_i//2))*_GCVS]]
    _StepList       += [[-_CrossThickness/2,-_CrossThickness/2]]
    return _StepList

def PathToPadFirstFractionSteps(_StartingPoint,_PadPoint1,_PadPoint2,_EndPoint,_Fraction):
    _Point1         = IntermediatePoint(_StartingPoint,_PadPoint1,_Fraction)
    _Point2         = IntermediatePoint(_EndPoint,_PadPoint2,_Fraction)
    return PositionsToSteps([_StartingPoint,_Point1,_Point2,_EndPoint])

def PathToPadSecondFractionSteps(_StartingPoint,_PadPoint1,_PadPoint2,_EndPoint,_Fraction):
    _Point1         = IntermediatePoint(_StartingPoint,_PadPoint1,_Fraction)
    _Point2         = IntermediatePoint(_EndPoint,_PadPoint2,_Fraction)
    return PositionsToSteps([_StartingPoint,_Point1,_PadPoint1,_PadPoint2,_Point2,_Point1,_StartingPoint])

# File operations
def WriteStepsToFile(_StartingPoint,_StepList):
    global PointCounter
    _PointList      = [_StartingPoint]
    for _Step in _StepList:
        _PointList.append(PosSum(_PointList[len(_PointList)-1],_Step))
    for _Point in _PointList:
        File.write("%d	%d\n" % (_Point[0],_Point[1]))
        PointCounter+=1
    File.write("\n") 

#%% Derived parameters
TLMTopLeft              = CoordinateToPosition(TLMTopLeftCoordinate)
GreekTopLeft            = CoordinateToPosition(GreekTopLeftCoordinate)
HallTopLeft             = CoordinateToPosition(HallTopLeftCoordinate)


#%% Mask 2: Si/Ge/Al, etch down Al and Ge. DEVICES
PolygonFileName = 'Polygon_%d%02d%02d_%02dh%02d_%s_Mask2.txt' % (now.year, now.month, now.day, now.hour, now.minute, DeviceType)
File            = open(PolygonFileName,"w+")
PointCounter    = 0

# Empty  list of points to contact the pads to (ends of the arms sticking out bottom and top)
ArmContactList  = [[[[[0,0],[0,0]] for ArmIterator in range(0,8)] for HIterator in range(0,HorizontalCells)] for VIterator in range(0,VerticalCells)]

# --- TLM --- #
for MesaIterator in range(0,len(WMesaList)):
    for ContactIterator in range(0,len(LContactList)):
        TotalLength     = 2*LEnd+8*LContactList[ContactIterator]+sum(LChannelList)
        DeviceCenter    = PosListSum([TLMTopLeft,StepFactor(ContactSpacing,ContactIterator),StepFactor(MesaSpacing,MesaIterator)])
        ReferencePoint  = PosSum(DeviceCenter,[-TotalLength/2,WMesaList[MesaIterator]/2])

        # Edge that sticks out on the left
        StepList        = [[0,-WMesaList[MesaIterator]],
                            [LEnd,0],
                            [0,-LArm]]

        # Calculation of e-beam part of trace to the pads
        StartOfTrace    = StepsToPosition(ReferencePoint,StepList)
        PadPosition1    = PosSum(DeviceCenter,PadPositions[0][0])
        PadPosition2    = PosSum(DeviceCenter,PadPositions[0][1])
        EndOfTrace      = [StartOfTrace[0]+LContactList[ContactIterator],StartOfTrace[1]-LContactList[ContactIterator]]
        ArmContactList[MesaIterator][ContactIterator][0]=[StartOfTrace,EndOfTrace]

        StepList        += PathToPadFirstFractionSteps(StartOfTrace,PadPosition1,PadPosition2,EndOfTrace,Litho1Fraction)
        StepList        += [[0,LArm+LContactList[ContactIterator]]]

        # Arms below
        for ChannelIterator in range(0,len(LChannelList),2):
            if len(LChannelList)%2 or ChannelIterator<len(LChannelList)-1:
                StepList    += [[LChannelList[ChannelIterator],0],
                                 [LContactList[ContactIterator],0]]
            if ChannelIterator<len(LChannelList)-2:
                StepList.append([LChannelList[ChannelIterator+1],0])
                StepList    += [[0,-(1+(ChannelIterator==0 or ChannelIterator==2))*LArm]] # Long arms in the middle

                # Calculation of e-beam part of trace to the pads
                StartOfTrace= StepsToPosition(ReferencePoint,StepList)
                PadPosition1= PosSum(DeviceCenter,PadPositions[1+ChannelIterator/2][0])
                PadPosition2= PosSum(DeviceCenter,PadPositions[1+ChannelIterator/2][1])
                EndOfTrace  = [StartOfTrace[0]+LContactList[ContactIterator],StartOfTrace[1]-(1-ChannelIterator/2)*LContactList[ContactIterator]]
                ArmContactList[MesaIterator][ContactIterator][1+ChannelIterator/2]=[StartOfTrace,EndOfTrace]
        
                StepList        += PathToPadFirstFractionSteps(StartOfTrace,PadPosition1,PadPosition2,EndOfTrace,Litho1Fraction)
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

            # Calculation of e-beam part of trace to the pads
            StartOfTrace    = StepsToPosition(ReferencePoint,StepList)
            PadPosition1    = PosSum(DeviceCenter,PadPositions[7-ChannelIterator/2][0])
            PadPosition2    = PosSum(DeviceCenter,PadPositions[7-ChannelIterator/2][1])
            EndOfTrace      = [StartOfTrace[0]-LContactList[ContactIterator],StartOfTrace[1]+(ChannelIterator/2-1-(ChannelIterator==6))*LContactList[ContactIterator]]
            ArmContactList[MesaIterator][ContactIterator][7-ChannelIterator/2]=[StartOfTrace,EndOfTrace]
    
            StepList        += PathToPadFirstFractionSteps(StartOfTrace,PadPosition1,PadPosition2,EndOfTrace,Litho1Fraction)
            StepList        += [[0,-(1+(ChannelIterator==2 or ChannelIterator==4))*LArm-(ChannelIterator/2-1-(ChannelIterator==6))*LContactList[ContactIterator]],
                                 [-LChannelList[ChannelIterator],0],
                                 [-LContactList[ContactIterator],0]]
            if ChannelIterator>1:
                StepList.append([-LChannelList[ChannelIterator-1],0])
        StepList.append([-LEnd,0])

        # Go back to the origin
        StepList        += [PosDiff(StepsToPosition(ReferencePoint,StepList),[0,0])]
        # Convert list of steps to list of points and write them to the file
        WriteStepsToFile(ReferencePoint,StepList)

# --- Hall bars --- #
for MesaIterator in range(0,len(HallWMesaList)):
    for LayerIterator in range(0,len(HallLayerList)):
        TotalLength     = HallLMesa
        DeviceCenter    = PosListSum([HallTopLeft,StepFactor(HallLayerSpacing,LayerIterator),StepFactor(HallMesaSpacing,MesaIterator)])
        ReferencePoint  = PosSum(DeviceCenter,[-TotalLength/2,HallWMesaList[MesaIterator]/2])
        HallLEnd        = (HallLMesa-2*HallArmSpacing-HallWContact)/2
        
        StepList        = [[0,0]]
        for SideIterator in range(0,2):   # Bottom and top are symmetric, just (-1)**SideIterator times all steps
            StartOfTrace= StepsToPosition(ReferencePoint,StepList)
            PadPosition1= PosSum(DeviceCenter,PadPositions[4*SideIterator][0])
            PadPosition2= PosSum(DeviceCenter,PadPositions[4*SideIterator][1])
            EndOfTrace  = [StartOfTrace[0],StartOfTrace[1]-(-1)**SideIterator*HallWMesaList[MesaIterator]]
            ArmContactList[MesaIterator][len(LContactList)+LayerIterator][4*SideIterator]=[StartOfTrace,EndOfTrace]
            StepList    += PathToPadFirstFractionSteps(StartOfTrace,PadPosition1,PadPosition2,EndOfTrace,Litho1Fraction)
            StepList    += [[(-1)**SideIterator*HallLEnd,0]]                
            for ArmIterator in range(0,3):
                StepList    += [[0,-(-1)**SideIterator*(HallLArm+0.5*ArmIterator*HallWContact)]]

                # Calculation of e-beam part of trace to the pads
                StartOfTrace= StepsToPosition(ReferencePoint,StepList)
                PadPosition1= PosSum(DeviceCenter,PadPositions[4*SideIterator+ArmIterator+1][0])
                PadPosition2= PosSum(DeviceCenter,PadPositions[4*SideIterator+ArmIterator+1][1])
                EndOfTrace  = [StartOfTrace[0]+(-1)**SideIterator*HallWContact,StartOfTrace[1]-(-1)**SideIterator*(1-ArmIterator)*HallWContact]
                ArmContactList[MesaIterator][len(LContactList)+LayerIterator][4*SideIterator+ArmIterator+1]=[StartOfTrace,EndOfTrace]
                StepList    += PathToPadFirstFractionSteps(StartOfTrace,PadPosition1,PadPosition2,EndOfTrace,Litho1Fraction)
                StepList    += [[0,(-1)**SideIterator*(HallLArm-0.5*(ArmIterator-2)*HallWContact)],[(-1)**SideIterator*(ArmIterator<2)*(HallArmSpacing-HallWContact),0]]
            StepList        += [[(-1)**SideIterator*HallLEnd,0]]

        # Go back to the origin
        StepList            += [PosDiff(StepsToPosition(ReferencePoint,StepList),[0,0])]
        # Convert list of steps to list of points and write them to the file
        WriteStepsToFile(ReferencePoint,StepList)

# --- Greek crosses --- #
for WidthIterator in range(0,len(GreekWidthList)):
    for LengthIterator in range(0,len(GreekLengthList)):
        for LayerIterator in range(0,2):
            TotalLength     = GreekLengthList[LengthIterator]
            DeviceCenter    = PosListSum([GreekTopLeft,StepFactor(GreekLengthSpacing,LengthIterator),StepFactor(GreekWidthSpacing,WidthIterator),[0,(1-2*LayerIterator)*PadSpacing/3]])
            ReferencePoint  = DeviceCenter
            
            GreekContactNrs = [[6,7,0,5],[4,1,2,3]]
            Thickness       = GreekWidthList[WidthIterator]
            CS              = GreekLengthList[LengthIterator]-Thickness/2
            StepList        = [[Thickness/2,Thickness/2]]
            for _i in range(0,4,1):
                StepList    += [[(_i%2)*((-1)**(1+_i//2))*CS,(1-_i%2)*((-1)**(_i//2))*CS]]
                StartOfTrace= StepsToPosition(ReferencePoint,StepList)
                PadPosition1= PosListSum([DeviceCenter,[0,(2*LayerIterator-1)*PadSpacing/3],PadPositions[GreekContactNrs[LayerIterator][_i]][0]])
                PadPosition2= PosListSum([DeviceCenter,[0,(2*LayerIterator-1)*PadSpacing/3],PadPositions[GreekContactNrs[LayerIterator][_i]][1]])
                EndOfTrace  = [StartOfTrace[0]-(1-_i%2)*(-1)**(_i//2)*Thickness,StartOfTrace[1]-((_i%2)*(-1)**(_i//2)+(LayerIterator==0 and _i==2)-(LayerIterator==1 and _i==0))*Thickness]
                ArmContactList[GreekTopLeftCoordinate[1]+WidthIterator][GreekTopLeftCoordinate[0]+LengthIterator][GreekContactNrs[LayerIterator][_i]]=[StartOfTrace,EndOfTrace]
                StepList    += PathToPadFirstFractionSteps(StartOfTrace,PadPosition1,PadPosition2,EndOfTrace,Litho1Fraction)
                StepList    += [[0,((LayerIterator==0 and _i==2)-(LayerIterator==1 and _i==0))*Thickness],[(_i%2)*((-1)**(_i//2))*CS,(1-_i%2)*((-1)**(1+_i//2))*CS]]
            StepList        += [[-Thickness/2,-Thickness/2]]

            # Go back to the origin
            StepList            += [PosDiff(StepsToPosition(ReferencePoint,StepList),[0,0])]
            # Convert list of steps to list of points and write them to the file
            WriteStepsToFile(ReferencePoint,StepList)

File.close()
print('File:\n    '+PolygonFileName+'\n')
print("Total number of points:\n    %d" % PointCounter)

#%% Mask 1: Si/Ge/Al, etch down Al and Ge. ALIGNMENT AND CONTACT PADS
PolygonFileName = 'Polygon_%d%02d%02d_%02dh%02d_%s_Mask1_ALIGNMENT_AND_CONTACT_PADS.txt' % (now.year, now.month, now.day, now.hour, now.minute, DeviceType)
File            = open(PolygonFileName,"w+")
PointCounter    = 0

# Alignment marks
StepList                = [[0,0]]
for Point in GlobalCrossPoints:
    StepList            += MoveTo([0,0],StepList,Point)+CrossSteps(GlobalCrossThickness,GlobalCrossWidth)
    StepList            += MoveTo([0,0],StepList,[0,0])
for Point in GlobalLongCrossPoints:
    StepList            += MoveTo([0,0],StepList,Point)+LongCrossSteps(GlobalCrossThickness,GlobalCrossWidth,GlobalCrossLongLength)
    StepList            += MoveTo([0,0],StepList,[0,0])
for GroupPoint in ChipMarkGroupPoints:
    for SetPoint in ChipMarkSetPoints:
        for Point in ChipMarkPoints:
            StepList        += MoveTo([0,0],StepList,PosListSum([GroupPoint,SetPoint,Point]))+CrossSteps(ChipCrossThickness,ChipCrossWidth)
            StepList        += MoveTo([0,0],StepList,[0,0])
WriteStepsToFile([0,0],StepList)

# Contact pads
for VIterator in range(0,VerticalCells):
    for HIterator in range(0,HorizontalCells):
        DeviceCenter    = PosListSum([TLMTopLeft,StepFactor(ContactSpacing,HIterator),StepFactor(MesaSpacing,VIterator)])
        ReferencePoint  = DeviceCenter
        StepList        = []

        # Contact pads
        for PadCenter in ContactPadList:
            StepList    += PositionsToSteps([ReferencePoint,
                            PosSum(PosSum(DeviceCenter,PadCenter),[PadWidth/2,PadWidth/2]),
                            PosSum(PosSum(DeviceCenter,PadCenter),[-PadWidth/2,PadWidth/2]),
                            PosSum(PosSum(DeviceCenter,PadCenter),[-PadWidth/2,-PadWidth/2]),
                            PosSum(PosSum(DeviceCenter,PadCenter),[PadWidth/2,-PadWidth/2]),
                            PosSum(PosSum(DeviceCenter,PadCenter),[PadWidth/2,PadWidth/2]),
                            ReferencePoint])
        for PadIterator in range(0,8):
            StartOfTrace= ArmContactList[VIterator][HIterator][PadIterator][0]
            EndOfTrace  = ArmContactList[VIterator][HIterator][PadIterator][1]
            PadPoint1   = PosSum(DeviceCenter,PadPositions[PadIterator][0])
            PadPoint2   = PosSum(DeviceCenter,PadPositions[PadIterator][1])
            
            StepList    += PositionsToSteps([ReferencePoint,StartOfTrace])
            StepList    += PathToPadSecondFractionSteps(StartOfTrace,PadPoint1,PadPoint2,EndOfTrace,Litho2Fraction)
            StepList    += PositionsToSteps([StartOfTrace,ReferencePoint])
        
        # Go back to the origin
        StepList        += [PosDiff(StepsToPosition(ReferencePoint,StepList),[0,0])]
        # Convert list of steps to list of points and write them to the file
        WriteStepsToFile(ReferencePoint,StepList)

File.close()
print('File:\n    '+PolygonFileName+'\n')
print("Total number of points:\n    %d" % PointCounter)

#%% Mask 3: Si/Ge/Al, etch down only Al.
PolygonFileName = 'Polygon_%d%02d%02d_%02dh%02d_%s_Mask3.txt' % (now.year, now.month, now.day, now.hour, now.minute, DeviceType)
File            = open(PolygonFileName,"w+")
PointCounter    = 0

# TLM
for MesaIterator in range(0,len(WMesaList)):
    for ContactIterator in range(0,len(LContactList)):
        TotalLength     = 2*LEnd+2*TLMMargin+8*LContactList[ContactIterator]+sum(LChannelList)
        DeviceCenter    = PosListSum([TLMTopLeft,StepFactor(ContactSpacing,ContactIterator),StepFactor(MesaSpacing,MesaIterator)])
        ReferencePoint  = PosSum(DeviceCenter,[-TotalLength/2,WMesaList[MesaIterator]/2+TLMMargin])
        
        # Edge that sticks out on the left
        StepList        = [[0,-WMesaList[MesaIterator]-2*TLMMargin],
                            [LEnd+TLMMargin,0],
                            [0,WMesaList[MesaIterator]+2*TLMMargin],
                            [-LEnd-TLMMargin,0],
                            [LEnd+TLMMargin+LContactList[ContactIterator],0]]

        # Arms below
        for ChannelIterator in range(0,len(LChannelList)):
            StepList    += [[0,-WMesaList[MesaIterator]-2*TLMMargin],
                            [LChannelList[ChannelIterator],0],
                            [0,WMesaList[MesaIterator]+2*TLMMargin],
                            [-LChannelList[ChannelIterator],0],
                            [LChannelList[ChannelIterator]+LContactList[ContactIterator],0]]

        # Edge that sticks out on the right
        StepList        += [[0,-WMesaList[MesaIterator]-2*TLMMargin],
                            [LEnd+TLMMargin,0],
                            [0,WMesaList[MesaIterator]+2*TLMMargin],
                            [-LEnd-TLMMargin,0],
                            [LEnd+TLMMargin-TotalLength,0]]

        # Go back to the origin
        StepList        += [PosDiff(StepsToPosition(ReferencePoint,StepList),[0,0])]
        # Convert list of steps to list of points and write them to the file
        WriteStepsToFile(ReferencePoint,StepList)

# Hall bars
for MesaIterator in range(0,len(HallWMesaList)):
    for LayerIterator in range(0,len(HallLayerList)):
        if HallLayerList[LayerIterator]:
            TotalLength     = HallLMesa+2*HallMargin
            DeviceCenter    = PosListSum([HallTopLeft,StepFactor(HallLayerSpacing,LayerIterator),StepFactor(HallMesaSpacing,MesaIterator)])
            ReferencePoint  = PosSum(DeviceCenter,[-TotalLength/2,HallWMesaList[MesaIterator]/2+HallLArm])
            StepList        = []
            for x in [0,1]:
                StepList    += StepListFactor([[0,-HallWMesaList[MesaIterator]-2*HallMargin-2*HallLArm],[TotalLength,0]],
                                (-1)**x)
    
            # Go back to the origin
            StepList        += [PosDiff(StepsToPosition(ReferencePoint,StepList),[0,0])]
            # Convert list of steps to list of points and write them to the file
            WriteStepsToFile(ReferencePoint,StepList)

# Greek crosses
for WidthIterator in range(0,len(GreekWidthList)):
    for LengthIterator in range(0,len(GreekLengthList)):
        for LayerIterator in range(0,2):
            if GreekLayerList[LayerIterator]:
                TotalLength     = 2*GreekLengthList[LengthIterator]+2*GreekMargin
                DeviceCenter    = PosListSum([GreekTopLeft,StepFactor(GreekLengthSpacing,LengthIterator),StepFactor(GreekWidthSpacing,WidthIterator),[0,(1-2*LayerIterator)*PadSpacing/3]])
                ReferencePoint  = PosSum(DeviceCenter,[-TotalLength/2,TotalLength/2])
                StepList        = []
                for x in [0,1]:
                    StepList    += StepListFactor([[0,-TotalLength],[TotalLength,0]],
                                    (-1)**x)
        
                # Go back to the origin
                StepList        += [PosDiff(StepsToPosition(ReferencePoint,StepList),[0,0])]
                # Convert list of steps to list of points and write them to the file
                WriteStepsToFile(ReferencePoint,StepList)  
  
File.close()
print('File:\n    '+PolygonFileName+'\n')
print("Total number of points:\n    %d" % PointCounter)