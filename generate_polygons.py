# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 08:59:33 2019

@author: TV255140
"""

from math import *
import numpy as np
from datetime import datetime
now = datetime.now()

DeviceType              = "All"


#%% Chip dimensions

HorizontalCells         = 6
HorizontalCellSpacing   = 1000000
VerticalCells           = 8
VerticalCellSpacing     = -1000000
ChipWidth               = HorizontalCells*1000000
ChipLength              = VerticalCells*1000000

# Alignment marks
GlobalCrossThickness    = 10000
GlobalCrossWidth        = 500000
GlobalCrossLongLength   = 2000000
GlobalCrossPoints       = [[-ChipWidth/2.,ChipLength/2.],                             [ChipWidth/2.,ChipLength/2.],
                           [-ChipWidth/2.,0],                                        [ChipWidth/2.,0],
                           [-ChipWidth/2.,-ChipLength/2.],    [0,-ChipLength/2.],      [ChipWidth/2.,-ChipLength/2.],]
GlobalLongCrossPoints   = [[0,ChipLength/2.]]
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
ChipMarkPoints          = [[-ChipMarkSpacing/2.,ChipMarkSpacing/2.],  [ChipMarkSpacing/2.,ChipMarkSpacing/2.],
                           [-ChipMarkSpacing/2.,-ChipMarkSpacing/2.], [ChipMarkSpacing/2.,-ChipMarkSpacing/2.]]

# Contact pads in a 3x3 grid
PadSpacing              = 300000
PadWidth                = 150000
WArm                    = 50000                         # Width of arm at the contact
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


#%% Device dimensions

# --- TLM --- #
TLMHorizontalCoordinates        = [0,2,3,5]
TLMVerticalCoordinates          = [0,2,4,6]
WMesaList                       = [1000,2000,5000,10000]        # All lengths in nanometers
LContactList                    = [50,100,200,500]              # 4 contact lengths per chip
LChannelList                    = [50,100,150,200,500,1000,2000]# 8 contacts on one TLM bar, 7 steps
LEnd                            = 1000                          # Length that edges of Si/Ge stick out left and right
LArm                            = 750                           # Vertical length that arms stick out before bending to contacts
TLMMargin                       = 100                           # Spill of mask 2 around the mesa

# --- Hall bars --- #
HallHorizontalCoordinates       = [1,4]
HallVerticalCoordinates         = [0,2,4,6]
HallWMesaList                   = [5000,10000,20000,50000]
HallLayerList                   = [1,0]
HallArmSpacing                  = 50000
HallLMesa                       = 200000
HallWContact                    = 2000
HallLArm                        = 3000
HallMargin                      = 0

# --- Greek crosses --- #
GreekHorizontalCoordinates      = [1,4]
GreekVerticalCoordinates        = [1,3,5,7]
GreekWidthList                  = [2000,5000,10000,20000]
GreekLengthList                 = [20000,50000]
GreekLayerList                  = [1,0]
GreekMargin                     = 250

# --- Meanders --- #
MeanderHorizontalCoordinates    = [0,2,3,5]
MeanderVerticalCoordinates      = [1,3,5,7]
MeanderWidthList                = [50,100,150,200]
MeanderLengthList               = [500000,500000,1000000,1000000]
MeanderLayerList                = [0,1,0,1]     # The grid has rows in width and columns in length, layer is applied per column
MeanderLayerMargin              = 250
MeanderContactWidth             = 50
MeanderContactLength            = 150
MeanderTurnRadius               = 25000
MeanderNDoubleTurns             = 2
MeanderNTurnPoints              = 64


#%% Parameter calculations
def CoordinateToPosition(Coordinate):
    return [((0.5-HorizontalCells/2.)+Coordinate[0])*HorizontalCellSpacing,((Coordinate[1]-VerticalCells/2.+0.5)*VerticalCellSpacing)]

# Empty  list of points to contact the pads to (ends of the arms sticking out bottom and top)
ArmContactList  = [[[[CoordinateToPosition([HIterator,VIterator]),
                      CoordinateToPosition([HIterator,VIterator])] 
                      for ArmIterator in range(0,8)] for HIterator in range(0,HorizontalCells)] 
                      for VIterator in range(0,VerticalCells)]

#%% Function definitions

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
    _CS             = _CrossWidth-_CrossThickness/2.            # Cross Step size
    _StepList       = [[_CrossThickness/2.,_CrossThickness/2.]]
    for _i in range(0,4,1):
        _StepList   += [[(_i%2)*((-1)**(1+_i//2.))*_CS,(1-_i%2)*((-1)**(_i//2.))*_CS],
                        [(1-_i%2)*((-1)**(1+_i//2.))*_CrossThickness,(_i%2)*((-1)**(1+_i//2.))*_CrossThickness],
                        [(_i%2)*((-1)**(_i//2.))*_CS,(1-_i%2)*((-1)**(1+_i//2.))*_CS]]
    _StepList       += [[-_CrossThickness/2.,-_CrossThickness/2.]]
    return _StepList

def LongCrossSteps(_CrossThickness,_CrossWidth,_CrossLength):
    _GCHS           = _CrossLength-_CrossThickness            # Global Cross Horizontal Step size
    _GCVS           = _CrossWidth-_CrossThickness
    _StepList       = [[_CrossThickness/2.,_CrossThickness/2.]]
    for _i in range(0,4,1):
        _StepList   += [[(_i%2)*((-1)**(1+_i//2.))*_GCHS,(1-_i%2)*((-1)**(_i//2.))*_GCVS],
                        [(1-_i%2)*((-1)**(1+_i//2.))*_CrossThickness,(_i%2)*((-1)**(1+_i//2.))*_CrossThickness],
                        [(_i%2)*((-1)**(_i//2.))*_GCHS,(1-_i%2)*((-1)**(1+_i//2.))*_GCVS]]
    _StepList       += [[-_CrossThickness/2.,-_CrossThickness/2.]]
    return _StepList

def PathToPadFirstFraction(_PadNr,_StartingPoint,_EndPoint):
    global HorizontalIterator, VerticalIterator, ArmContactList
    _DeviceCenter   = CoordinateToPosition([HorizontalIterator,VerticalIterator])
    _PadPoint1      = PosSum(_DeviceCenter,PadPositions[int(_PadNr)][0])
    _PadPoint2      = PosSum(_DeviceCenter,PadPositions[int(_PadNr)][1])
    _Point1         = IntermediatePoint(_StartingPoint,_PadPoint1,Litho1Fraction)
    _Point2         = IntermediatePoint(_EndPoint,_PadPoint2,Litho1Fraction)
    ArmContactList[VerticalIterator][HorizontalIterator][int(_PadNr)]=[_StartingPoint,_EndPoint]
    return PositionsToSteps([_StartingPoint,_Point1,_Point2,_EndPoint]) 

def PathToPadSecondFraction(_PadNr):
    global HorizontalIterator, VerticalIterator, ArmContactList
    _StartingPoint  = ArmContactList[VerticalIterator][HorizontalIterator][int(_PadNr)][0]
    _EndPoint       = ArmContactList[VerticalIterator][HorizontalIterator][int(_PadNr)][1]
    _DeviceCenter   = CoordinateToPosition([HorizontalIterator,VerticalIterator])
    _PadPoint1      = PosSum(_DeviceCenter,PadPositions[int(_PadNr)][0])
    _PadPoint2      = PosSum(_DeviceCenter,PadPositions[int(_PadNr)][1])
    _Point1         = IntermediatePoint(_StartingPoint,_PadPoint1,Litho2Fraction)
    _Point2         = IntermediatePoint(_EndPoint,_PadPoint2,Litho2Fraction)
    return PositionsToSteps([_DeviceCenter,_StartingPoint,_Point1,_PadPoint1,_PadPoint2,_Point2,_Point1,_StartingPoint,_DeviceCenter])

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

#%% Functions to draw structures

# --- TLM --- #
def DrawTLM():
    global HorizontalIterator, VerticalIterator, ArmContactList
    for MesaIterator in range(0,len(_WMesaList)):
        for ContactIterator in range(0,len(_LContactList)):
            HorizontalIterator  = TLMTopLeftCoordinate[0]+ContactIterator
            VerticalIterator    = TLMTopLeftCoordinate[1]+MesaIterator
            TotalLength         = 2*LEnd+8*_LContactList[ContactIterator]+sum(LChannelList)
            DeviceCenter        = CoordinateToPosition([HorizontalIterator,VerticalIterator])
            ReferencePoint      = PosSum(DeviceCenter,[-TotalLength/2.,_WMesaList[MesaIterator]/2.])
    
            # Edge that sticks out on the left
            StepList        = [[0,-_WMesaList[MesaIterator]],
                                [LEnd,0],
                                [0,-LArm]]
    
            # Calculation of first part of trace to the pads
            StartOfTrace    = StepsToPosition(ReferencePoint,StepList)
            EndOfTrace      = [StartOfTrace[0]+_LContactList[ContactIterator],StartOfTrace[1]-_LContactList[ContactIterator]]
            StepList        += PathToPadFirstFraction(0,StartOfTrace,EndOfTrace)
            StepList        += [[0,LArm+_LContactList[ContactIterator]]]
    
            # Arms below
            for ChannelIterator in range(0,len(LChannelList),2):
                if len(LChannelList)%2 or ChannelIterator<len(LChannelList)-1:
                    StepList    += [[LChannelList[ChannelIterator],0],
                                     [_LContactList[ContactIterator],0]]
                if ChannelIterator<len(LChannelList)-2:
                    StepList.append([LChannelList[ChannelIterator+1],0])
                    StepList    += [[0,-(1+(ChannelIterator==0 or ChannelIterator==2))*LArm]] # Long arms in the middle
    
                    # Calculation of e-beam part of trace to the pads
                    StartOfTrace= StepsToPosition(ReferencePoint,StepList)
                    EndOfTrace  = [StartOfTrace[0]+_LContactList[ContactIterator],StartOfTrace[1]-(1-ChannelIterator/2.)*_LContactList[ContactIterator]]
                    StepList    += PathToPadFirstFraction(1+ChannelIterator/2,StartOfTrace,EndOfTrace)
                    StepList    += [[0,(1+(ChannelIterator==0 or ChannelIterator==2))*LArm+(1-ChannelIterator/2.)*_LContactList[ContactIterator]]]  # Long arms in the middle
    
            # Edge that sticks out on the right
            StepList            += [[LEnd,0],
                                    [0,_WMesaList[MesaIterator]],
                                    [-LEnd,0]]
            if (len(LChannelList)+1)%2:     # If there is an odd number of arms, take a step first
                StepList.append([-LChannelList[len(LChannelList)-1]-_LContactList[ContactIterator],0])
    
            # Arms above
            for ChannelIterator in range(len(LChannelList)-len(LChannelList)%2,-1,-2):
                StepList        += [[0,(1+(ChannelIterator==2 or ChannelIterator==4))*LArm]]    # Long arms in the middle
    
                # Calculation of e-beam part of trace to the pads
                StartOfTrace    = StepsToPosition(ReferencePoint,StepList)
                EndOfTrace      = [StartOfTrace[0]-_LContactList[ContactIterator],StartOfTrace[1]+(ChannelIterator/2.-1-(ChannelIterator==6))*_LContactList[ContactIterator]]
                StepList        += PathToPadFirstFraction(7-ChannelIterator/2,StartOfTrace,EndOfTrace)
                StepList        += [[0,-(1+(ChannelIterator==2 or ChannelIterator==4))*LArm-(ChannelIterator/2.-1-(ChannelIterator==6))*_LContactList[ContactIterator]],
                                     [-LChannelList[ChannelIterator],0],
                                     [-_LContactList[ContactIterator],0]]
                if ChannelIterator>1:
                    StepList.append([-LChannelList[ChannelIterator-1],0])
            StepList.append([-LEnd,0])
    
            # Go back to the origin
            StepList        += [PosDiff(StepsToPosition(ReferencePoint,StepList),[0,0])]
            # Convert list of steps to list of points and write them to the file
            WriteStepsToFile(ReferencePoint,StepList)

# --- Hall bars --- #
def DrawHallBars():
    global HorizontalIterator, VerticalIterator, ArmContactList
    for MesaIterator in range(0,len(_HallWMesaList)):
        for LayerIterator in range(0,len(_HallLayerList)):
            TotalLength         = HallLMesa
            HorizontalIterator  = HallTopLeftCoordinate[0]+LayerIterator
            VerticalIterator    = HallTopLeftCoordinate[1]+MesaIterator
            DeviceCenter        = CoordinateToPosition([HorizontalIterator,VerticalIterator])
            ReferencePoint      = PosSum(DeviceCenter,[-TotalLength/2.,_HallWMesaList[MesaIterator]/2.])
            HallLEnd            = (HallLMesa-2*HallArmSpacing-HallWContact)/2.
            
            StepList            = [[0,0]]
            for SideIterator in range(0,2):   # Bottom and top are symmetric, just (-1)**SideIterator times all steps
                StartOfTrace    = StepsToPosition(ReferencePoint,StepList)
                EndOfTrace      = [StartOfTrace[0],StartOfTrace[1]-(-1)**SideIterator*_HallWMesaList[MesaIterator]]
                StepList        += PathToPadFirstFraction(4*SideIterator,StartOfTrace,EndOfTrace)
                StepList        += [[(-1)**SideIterator*HallLEnd,0]]                
                for ArmIterator in range(0,3):
                    StepList    += [[0,-(-1)**SideIterator*(HallLArm+0.5*ArmIterator*HallWContact)]]
    
                    # Calculation of e-beam part of trace to the pads
                    StartOfTrace= StepsToPosition(ReferencePoint,StepList)
                    EndOfTrace  = [StartOfTrace[0]+(-1)**SideIterator*HallWContact,StartOfTrace[1]-(-1)**SideIterator*(1-ArmIterator)*HallWContact]
                    StepList    += PathToPadFirstFraction(4*SideIterator+ArmIterator+1,StartOfTrace,EndOfTrace)
                    StepList    += [[0,(-1)**SideIterator*(HallLArm-0.5*(ArmIterator-2)*HallWContact)],[(-1)**SideIterator*(ArmIterator<2)*(HallArmSpacing-HallWContact),0]]
                StepList        += [[(-1)**SideIterator*HallLEnd,0]]
    
            # Go back to the origin
            StepList            += [PosDiff(StepsToPosition(ReferencePoint,StepList),[0,0])]
            # Convert list of steps to list of points and write them to the file
            WriteStepsToFile(ReferencePoint,StepList)

# --- Greek crosses --- #
def DrawGreekCrosses():
    global HorizontalIterator, VerticalIterator, ArmContactList
    for WidthIterator in range(0,len(_GreekWidthList)):
        for LengthIterator in range(0,len(_GreekLengthList)):
            for LayerIterator in range(0,len(_GreekLayerList)):
                HorizontalIterator  = GreekTopLeftCoordinate[0]+LengthIterator
                VerticalIterator    = GreekTopLeftCoordinate[1]+WidthIterator
                DeviceCenter        = PosSum(CoordinateToPosition([HorizontalIterator,VerticalIterator]),[0,(1-2*LayerIterator)*PadSpacing/3])
                ReferencePoint      = DeviceCenter
                
                GreekContactNrs     = [[6,7,0,5],[4,1,2,3]]
                Thickness           = _GreekWidthList[WidthIterator]
                CS                  = _GreekLengthList[LengthIterator]-Thickness/2.
                StepList            = [[Thickness/2.,Thickness/2.]]
                for _i in range(0,4,1):
                    StepList    += [[(_i%2)*((-1)**(1+_i//2.))*CS,(1-_i%2)*((-1)**(_i//2.))*CS]]
                    StartOfTrace= StepsToPosition(ReferencePoint,StepList)
                    EndOfTrace  = [StartOfTrace[0]-(1-_i%2)*(-1)**(_i//2.)*Thickness,StartOfTrace[1]-((_i%2)*(-1)**(_i//2.)+(LayerIterator==0 and _i==2)-(LayerIterator==1 and _i==0))*Thickness]
                    StepList    += PathToPadFirstFraction(GreekContactNrs[LayerIterator][_i],StartOfTrace,EndOfTrace)
                    StepList    += [[0,((LayerIterator==0 and _i==2)-(LayerIterator==1 and _i==0))*Thickness],[(_i%2)*((-1)**(_i//2.))*CS,(1-_i%2)*((-1)**(1+_i//2.))*CS]]
                StepList        += [[-Thickness/2.,-Thickness/2.]]
    
                # Go back to the origin
                StepList            += [PosDiff(StepsToPosition(ReferencePoint,StepList),[0,0])]
                # Convert list of steps to list of points and write them to the file
                WriteStepsToFile(ReferencePoint,StepList)

# --- Meanders --- #
def DrawMeanders():
    global HorizontalIterator, VerticalIterator, ArmContactList
    for WidthIterator in range(0,len(_MeanderWidthList)):
        for LengthIterator in range(0,len(_MeanderLengthList)):
            Length              = (_MeanderLengthList[LengthIterator]-(4*pi+2)*MeanderTurnRadius)/5
            Width               = _MeanderWidthList[WidthIterator]
            TotalLength         = Length+2*MeanderTurnRadius+Width
            TotalWidth          = 8*MeanderTurnRadius+Width
            HorizontalIterator  = MeanderTopLeftCoordinate[0]+LengthIterator
            VerticalIterator    = MeanderTopLeftCoordinate[1]+WidthIterator
            DeviceCenter        = CoordinateToPosition([HorizontalIterator,VerticalIterator])
            ReferencePoint      = PosSum(DeviceCenter,[TotalLength/2.-Width/2.,TotalWidth/2.-Width])
            
            StepList            = [[0,0]]
            
            # Move Length, then make a turn
            for Side in range(0,2):
    
                # First and last contacts
                StepList            += [[(-1)**Side*Width/2.,0]]
                StartOfTrace        = StepsToPosition(ReferencePoint,StepList)
                EndOfTrace          = [StartOfTrace[0]-(-1)**Side*Width/2.,StartOfTrace[1]+(-1)**Side*Width]
                StepList            += PathToPadFirstFraction(5-4*Side,StartOfTrace,EndOfTrace)
                
                StepList            += [[-(-1)**Side*(MeanderTurnRadius+Length/2-MeanderContactWidth/2),0],
                                        [0,(-1)**Side*MeanderContactLength]]
                StartOfTrace        = StepsToPosition(ReferencePoint,StepList)
                EndOfTrace          = [StartOfTrace[0]-(-1)**Side*MeanderContactWidth,StartOfTrace[1]]
                StepList            += PathToPadFirstFraction(6-Side*4,StartOfTrace,EndOfTrace)
                StepList            += [[0,-(-1)**Side*MeanderContactLength],[-(-1)**Side*(Length/2-MeanderContactWidth/2),0]]
                for IDoubleTurn in range(0,MeanderNDoubleTurns):
                    for Turn in range(0,2):
                        TurnCenter      = PosSum(StepsToPosition(ReferencePoint,StepList),[0,-(-1)**Side*MeanderTurnRadius])
                        Radius          = MeanderTurnRadius+(-1)**Turn*Width/2.
                        TurnPoints      = [PosSum(TurnCenter,[0,(-1)**Side*Radius])]
                        for TurnPoint in range(1,MeanderNTurnPoints+1):
                            Angle       = (-1)**Side*pi*(0.5+1.*TurnPoint/MeanderNTurnPoints)
                            TurnPoints  += [PosSum(TurnCenter,
                                                  [(-1)**(Side+Turn)*cos(Angle)*Radius,sin(Angle)*Radius])]
                        if 1-Turn:
                            StepList    += PositionsToSteps(TurnPoints[0:int(len(TurnPoints)/2)]
                                            + [PosSum(TurnPoints[int(len(TurnPoints)/2)],
                                                [0,(-1)**Side*MeanderContactWidth/2])])
                            StepList    += [[-(-1)**Side*MeanderContactLength,0]]
                            StartOfTrace= StepsToPosition(ReferencePoint,StepList)
                            EndOfTrace  = [StartOfTrace[0]-(-1)**Side*0.5*(1-IDoubleTurn)*MeanderContactWidth,StartOfTrace[1]-(-1)**Side*MeanderContactWidth]
                            PadNr       = (Side==0 and IDoubleTurn==0)*7+(Side==1 and IDoubleTurn==1)*4+(Side==1 and IDoubleTurn==0)*3
                            StepList    += PathToPadFirstFraction(PadNr,StartOfTrace,EndOfTrace)
                            StepList    +=[[(-1)**Side*((0.5-0.5*IDoubleTurn)*MeanderContactWidth+MeanderContactLength),0]]
                            StepList    += PositionsToSteps([PosSum(TurnPoints[int(len(TurnPoints)/2)],
                                                [0,-(-1)**Side*MeanderContactWidth/2])]
                                            + TurnPoints[int(len(TurnPoints)/2)+1:len(TurnPoints)+1])
                        else:
                            StepList    += PositionsToSteps(TurnPoints)
                        StepList        += [[(-1)**(Turn+Side)*Length,0]]
                StepList            += [[-(-1)**Side*MeanderTurnRadius,0]]
            
            # Go back to the origin
            StepList            += [PosDiff(StepsToPosition(ReferencePoint,StepList),[0,0])]
            # Convert list of steps to list of points and write them to the file
            WriteStepsToFile(ReferencePoint,StepList)

# Alignment marks
def DrawAlignmentMarks():
    global HorizontalIterator, VerticalIterator, ArmContactList
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
def DrawContactPads():
    global HorizontalIterator, VerticalIterator, ArmContactList
    for VerticalIterator in range(0,VerticalCells):
        for HorizontalIterator in range(0,HorizontalCells):
            DeviceCenter    = CoordinateToPosition([HorizontalIterator,VerticalIterator])
            ReferencePoint  = DeviceCenter
            StepList        = []
    
            # Contact pads
            for PadCenter in ContactPadList:
                StepList    += PositionsToSteps([ReferencePoint,
                                PosSum(PosSum(DeviceCenter,PadCenter),[PadWidth/2.,PadWidth/2.]),
                                PosSum(PosSum(DeviceCenter,PadCenter),[-PadWidth/2.,PadWidth/2.]),
                                PosSum(PosSum(DeviceCenter,PadCenter),[-PadWidth/2.,-PadWidth/2.]),
                                PosSum(PosSum(DeviceCenter,PadCenter),[PadWidth/2.,-PadWidth/2.]),
                                PosSum(PosSum(DeviceCenter,PadCenter),[PadWidth/2.,PadWidth/2.]),
                                ReferencePoint])
            for PadIterator in range(0,8):
                StepList    += PathToPadSecondFraction(PadIterator)
            
            # Go back to the origin
            StepList        += [PosDiff(StepsToPosition(ReferencePoint,StepList),[0,0])]
            # Convert list of steps to list of points and write them to the file
            WriteStepsToFile(ReferencePoint,StepList)

# TLM Aluminium
def DrawTLMAlu():
    global HorizontalIterator, VerticalIterator, ArmContactList
    for MesaIterator in range(0,len(_WMesaList)):
        for ContactIterator in range(0,len(_LContactList)):
            HorizontalIterator  = TLMTopLeftCoordinate[0]+ContactIterator
            VerticalIterator    = TLMTopLeftCoordinate[1]+MesaIterator
            TotalLength         = 2*LEnd+2*TLMMargin+8*_LContactList[ContactIterator]+sum(LChannelList)
            DeviceCenter        = CoordinateToPosition([HorizontalIterator,VerticalIterator])
            ReferencePoint      = PosSum(DeviceCenter,[-TotalLength/2.,_WMesaList[MesaIterator]/2.+TLMMargin])
            
            # Edge that sticks out on the left
            StepList        = [[0,-_WMesaList[MesaIterator]-2*TLMMargin],
                                [LEnd+TLMMargin,0],
                                [0,_WMesaList[MesaIterator]+2*TLMMargin],
                                [-LEnd-TLMMargin,0],
                                [LEnd+TLMMargin+_LContactList[ContactIterator],0]]
    
            # Arms below
            for ChannelIterator in range(0,len(LChannelList)):
                StepList    += [[0,-_WMesaList[MesaIterator]-2*TLMMargin],
                                [LChannelList[ChannelIterator],0],
                                [0,_WMesaList[MesaIterator]+2*TLMMargin],
                                [-LChannelList[ChannelIterator],0],
                                [LChannelList[ChannelIterator]+_LContactList[ContactIterator],0]]
    
            # Edge that sticks out on the right
            StepList        += [[0,-_WMesaList[MesaIterator]-2*TLMMargin],
                                [LEnd+TLMMargin,0],
                                [0,_WMesaList[MesaIterator]+2*TLMMargin],
                                [-LEnd-TLMMargin,0],
                                [LEnd+TLMMargin-TotalLength,0]]
    
            # Go back to the origin
            StepList        += [PosDiff(StepsToPosition(ReferencePoint,StepList),[0,0])]
            # Convert list of steps to list of points and write them to the file
            WriteStepsToFile(ReferencePoint,StepList)

# Hall bars
def DrawHallBarsAlu():
    global HorizontalIterator, VerticalIterator, ArmContactList
    for MesaIterator in range(0,len(_HallWMesaList)):
        for LayerIterator in range(0,len(_HallLayerList)):
            if _HallLayerList[LayerIterator]:
                TotalLength     = HallLMesa+2*HallMargin
                HorizontalIterator  = HallTopLeftCoordinate[0]+LayerIterator
                VerticalIterator    = HallTopLeftCoordinate[1]+MesaIterator
                DeviceCenter        = CoordinateToPosition([HorizontalIterator,VerticalIterator])
                ReferencePoint      = PosSum(DeviceCenter,[-TotalLength/2.,_HallWMesaList[MesaIterator]/2.+HallLArm])
                StepList            = []
                for x in [0,1]:
                    StepList        += StepListFactor([[0,-_HallWMesaList[MesaIterator]-2*HallMargin-2*HallLArm],[TotalLength,0]],
                                    (-1)**x)
        
                # Go back to the origin
                StepList            += [PosDiff(StepsToPosition(ReferencePoint,StepList),[0,0])]
                # Convert list of steps to list of points and write them to the file
                WriteStepsToFile(ReferencePoint,StepList)

# Greek crosses
def DrawGreekCrossesAlu():
    global HorizontalIterator, VerticalIterator, ArmContactList
    for WidthIterator in range(0,len(_GreekWidthList)):
        for LengthIterator in range(0,len(_GreekLengthList)):
            for LayerIterator in range(0,len(_GreekLayerList)):
                if _GreekLayerList[LayerIterator]:
                    TotalLength         = _GreekLengthList[LengthIterator]
                    TotalWidth          = _GreekWidthList[WidthIterator]+2*GreekMargin
                    HorizontalIterator  = GreekTopLeftCoordinate[0]+LengthIterator
                    VerticalIterator    = GreekTopLeftCoordinate[1]+WidthIterator
                    DeviceCenter        = PosSum(CoordinateToPosition([HorizontalIterator,VerticalIterator]),[0,(1-2*LayerIterator)*PadSpacing/3])
                    ReferencePoint      = DeviceCenter
                    StepList            = CrossSteps(TotalWidth,TotalLength)
    
                    # Go back to the origin
                    StepList            += [PosDiff(StepsToPosition(ReferencePoint,StepList),[0,0])]
                    # Convert list of steps to list of points and write them to the file
                    WriteStepsToFile(ReferencePoint,StepList)

# Meanders
def DrawMeandersAlu():
    global HorizontalIterator, VerticalIterator, ArmContactList
    for WidthIterator in range(0,len(_MeanderWidthList)):
        for LengthIterator in range(0,len(_MeanderLengthList)):
            if _MeanderLayerList[LengthIterator]:
                Length              = (_MeanderLengthList[LengthIterator]-(4*pi+2)*MeanderTurnRadius)/5
                Width               = _MeanderWidthList[WidthIterator]+2*MeanderLayerMargin
                TotalLength         = Length+2*MeanderTurnRadius+Width
                TotalWidth          = 8*MeanderTurnRadius+Width
                HorizontalIterator  = MeanderTopLeftCoordinate[0]+LengthIterator
                VerticalIterator    = MeanderTopLeftCoordinate[1]+WidthIterator
                DeviceCenter        = CoordinateToPosition([HorizontalIterator,VerticalIterator])
                ReferencePoint      = PosSum(DeviceCenter,[TotalLength/2.-Width/2.,TotalWidth/2.])
                
                StepList            = [[0,0]]
                
                # Move Length, then make a turn
                for Side in range(0,2):
                    
                    StepList        += [[-(-1)**Side*(Length+MeanderTurnRadius),0]]
                    for IDoubleTurn in range(0,MeanderNDoubleTurns):
                        for Turn in range(0,2):
                            TurnCenter      = PosSum(StepsToPosition(ReferencePoint,StepList),[0,-(-1)**Side*MeanderTurnRadius])
                            Radius          = MeanderTurnRadius+(-1)**Turn*Width/2.
                            TurnPoints      = [PosSum(TurnCenter,[0,(-1)**Side*Radius])]
                            for TurnPoint in range(1,MeanderNTurnPoints+1):
                                Angle       = (-1)**Side*pi*(0.5+1.*TurnPoint/MeanderNTurnPoints)
                                TurnPoints  += [PosSum(TurnCenter,
                                                      [(-1)**(Side+Turn)*cos(Angle)*Radius,sin(Angle)*Radius])]
                            StepList    += PositionsToSteps(TurnPoints)
                            StepList    += [[(-1)**(Turn+Side)*Length,0]]
                    StepList            += [[-(-1)**Side*MeanderTurnRadius,0],[0,-(-1)**Side*Width]]
                    
                # Go back to the origin
                StepList            += [PosDiff(StepsToPosition(ReferencePoint,StepList),[0,0])]
                # Convert list of steps to list of points and write them to the file
                WriteStepsToFile(ReferencePoint,StepList)

#%% Mask 3: Si/Ge/Al, etch down Al and Ge. DEVICES
PolygonFileName = 'Polygon_%d%02d%02d_%02dh%02d_%s_Mask3_Al_Ge.txt' % (now.year, now.month, now.day, now.hour, now.minute, DeviceType)
File            = open(PolygonFileName,"w+")
PointCounter    = 0


# --- TLM --- #
for TLMVerticalIterator in range(0,len(TLMVerticalCoordinates)):
    for TLMHorizontalIterator in range(0,len(TLMHorizontalCoordinates)):
        TLMTopLeftCoordinate    = [TLMHorizontalCoordinates[TLMHorizontalIterator],TLMVerticalCoordinates[TLMVerticalIterator]]
        _WMesaList               = [WMesaList[TLMVerticalIterator]]
        _LContactList            = [LContactList[TLMHorizontalIterator]]
        DrawTLM()

# --- Hall bars --- #
for HallVerticalIterator in range(0,len(HallVerticalCoordinates)):
    for HallHorizontalIterator in range(0,len(HallHorizontalCoordinates)):
        HallTopLeftCoordinate   = [HallHorizontalCoordinates[HallHorizontalIterator],HallVerticalCoordinates[HallVerticalIterator]]
        _HallWMesaList           = [HallWMesaList[HallVerticalIterator]]
        _HallLayerList           = [HallLayerList[HallHorizontalIterator]]
        DrawHallBars()

# --- Greek crosses --- #
for GreekVerticalIterator in range(0,len(GreekVerticalCoordinates)):
    for GreekHorizontalIterator in range(0,len(GreekHorizontalCoordinates)):
        GreekTopLeftCoordinate  = [GreekHorizontalCoordinates[GreekHorizontalIterator],GreekVerticalCoordinates[GreekVerticalIterator]]
        _GreekWidthList          = [GreekWidthList[GreekVerticalIterator]]
        _GreekLengthList         = [GreekLengthList[GreekHorizontalIterator]]
        _GreekLayerList          = GreekLayerList
        DrawGreekCrosses()

# --- Meanders --- #
for MeanderVerticalIterator in range(0,len(MeanderVerticalCoordinates)):
    for MeanderHorizontalIterator in range(0,len(MeanderHorizontalCoordinates)):
        MeanderTopLeftCoordinate  = [MeanderHorizontalCoordinates[MeanderHorizontalIterator],MeanderVerticalCoordinates[MeanderVerticalIterator]]
        _MeanderWidthList          = [MeanderWidthList[MeanderVerticalIterator]]
        _MeanderLengthList         = [MeanderLengthList[MeanderHorizontalIterator]]
        _MeanderLayerList          = [MeanderLayerList[MeanderHorizontalIterator]]
        DrawMeanders()


File.close()
print('File:\n    '+PolygonFileName+'\n')
print("Total number of points:\n    %d" % PointCounter)

#%% Mask 1: Si/Ge/Al, etch down Al and Ge. ALIGNMENT MARKS
PolygonFileName = 'Polygon_%d%02d%02d_%02dh%02d_%s_Mask1_ALIGNMENT_MARKS.txt' % (now.year, now.month, now.day, now.hour, now.minute, DeviceType)
File            = open(PolygonFileName,"w+")
PointCounter    = 0

# Alignment marks
DrawAlignmentMarks()

File.close()
print('File:\n    '+PolygonFileName+'\n')
print("Total number of points:\n    %d" % PointCounter)


#%% Mask 2: Si/Ge/Al, etch down Al and Ge. CONTACT PADS
PolygonFileName = 'Polygon_%d%02d%02d_%02dh%02d_%s_Mask2_CONTACT_PADS.txt' % (now.year, now.month, now.day, now.hour, now.minute, DeviceType)
File            = open(PolygonFileName,"w+")
PointCounter    = 0

# Contact pads
DrawContactPads()

File.close()
print('File:\n    '+PolygonFileName+'\n')
print("Total number of points:\n    %d" % PointCounter)

#%% Mask 4: Si/Ge/Al, etch down only Al.
PolygonFileName = 'Polygon_%d%02d%02d_%02dh%02d_%s_Mask4_Al.txt' % (now.year, now.month, now.day, now.hour, now.minute, DeviceType)
File            = open(PolygonFileName,"w+")
PointCounter    = 0

# TLM
for TLMVerticalIterator in range(0,len(TLMVerticalCoordinates)):
    for TLMHorizontalIterator in range(0,len(TLMHorizontalCoordinates)):
        TLMTopLeftCoordinate    = [TLMHorizontalCoordinates[TLMHorizontalIterator],TLMVerticalCoordinates[TLMVerticalIterator]]
        _WMesaList               = [WMesaList[TLMVerticalIterator]]
        _LContactList            = [LContactList[TLMHorizontalIterator]]
        DrawTLMAlu()

# Hall bars
for HallVerticalIterator in range(0,len(HallVerticalCoordinates)):
    for HallHorizontalIterator in range(0,len(HallHorizontalCoordinates)):
        HallTopLeftCoordinate   = [HallHorizontalCoordinates[HallHorizontalIterator],HallVerticalCoordinates[HallVerticalIterator]]
        _HallWMesaList          = [HallWMesaList[HallVerticalIterator]]
        _HallLayerList          = [HallLayerList[HallHorizontalIterator]]
        DrawHallBarsAlu()

# Greek crosses
for GreekVerticalIterator in range(0,len(GreekVerticalCoordinates)):
    for GreekHorizontalIterator in range(0,len(GreekHorizontalCoordinates)):
        GreekTopLeftCoordinate  = [GreekHorizontalCoordinates[GreekHorizontalIterator],GreekVerticalCoordinates[GreekVerticalIterator]]
        _GreekWidthList          = [GreekWidthList[GreekVerticalIterator]]
        _GreekLengthList         = [GreekLengthList[GreekHorizontalIterator]]
        _GreekLayerList          = GreekLayerList
        DrawGreekCrossesAlu()

# Meanders
for MeanderVerticalIterator in range(0,len(MeanderVerticalCoordinates)):
    for MeanderHorizontalIterator in range(0,len(MeanderHorizontalCoordinates)):
        MeanderTopLeftCoordinate  = [MeanderHorizontalCoordinates[MeanderHorizontalIterator],MeanderVerticalCoordinates[MeanderVerticalIterator]]
        _MeanderWidthList          = [MeanderWidthList[MeanderVerticalIterator]]
        _MeanderLengthList         = [MeanderLengthList[MeanderHorizontalIterator]]
        _MeanderLayerList          = [MeanderLayerList[MeanderHorizontalIterator]]
        DrawMeandersAlu()

File.close()
print('File:\n    '+PolygonFileName+'\n')
print("Total number of points:\n    %d" % PointCounter)


#%% Output mask layout as a table

DeviceGrid    = [["" for x in range(HorizontalCells)] for y in range(VerticalCells)]

for TLMVerticalIterator in range(0,len(TLMVerticalCoordinates)):
    for TLMHorizontalIterator in range(0,len(TLMHorizontalCoordinates)):
        DeviceGrid[TLMVerticalCoordinates[TLMVerticalIterator]][TLMHorizontalCoordinates[TLMHorizontalIterator]]   = "T"

# Hall bars
for HallVerticalIterator in range(0,len(HallVerticalCoordinates)):
    for HallHorizontalIterator in range(0,len(HallHorizontalCoordinates)):
        DeviceGrid[HallVerticalCoordinates[HallVerticalIterator]][HallHorizontalCoordinates[HallHorizontalIterator]]   = "H"

# Greek crosses
for GreekVerticalIterator in range(0,len(GreekVerticalCoordinates)):
    for GreekHorizontalIterator in range(0,len(GreekHorizontalCoordinates)):
        DeviceGrid[GreekVerticalCoordinates[GreekVerticalIterator]][GreekHorizontalCoordinates[GreekHorizontalIterator]]   = "G"

# Meanders
for MeanderVerticalIterator in range(0,len(MeanderVerticalCoordinates)):
    for MeanderHorizontalIterator in range(0,len(MeanderHorizontalCoordinates)):
        DeviceGrid[MeanderVerticalCoordinates[MeanderVerticalIterator]][MeanderHorizontalCoordinates[MeanderHorizontalIterator]]   = "M"
        
for Device in DeviceGrid:
    print(", ".join([str(l).rjust(2) for l in Device]))