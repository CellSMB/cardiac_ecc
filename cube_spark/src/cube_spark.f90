!> \file
!> $Id: cube_spark.f90 1528 2010-09-21 01:32:29Z  $
!> \author Chris Bradley
!> \brief This is test example of calcium spark kinetics and diffusion in cardiac cells.
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>
!> The Original Code is OpenCMFE_
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand and University of Oxford, Oxford, United
!> Kingdom. Portions created by the University of Auckland and University
!> of Oxford are Copyright (C) 2007 by the University of Auckland and
!> the University of Oxford. All Rights Reserved.
!>
!> Contributor(s):
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>

!> \example cube_spark.f90
!! Example program which sets up a single calcium release site inside a cubic cardiac cell geometry .
!! \par Latest Builds:
!<

!> Main program
PROGRAM CUBE_SPARK

  USE OPENCMISS
  USE OpenCMISS_Iron
#ifndef NOMPIMOD
  USE MPI
#endif

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

  ! program parameters

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=3
  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=4
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=6


  INTEGER(CMISSIntg), PARAMETER :: CaEquationsSetUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: CaMaterialsFieldUserNumber=8
  INTEGER(CMISSIntg), PARAMETER :: CaFieldUserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: CaEquationsSetFieldUserNumber=10
  INTEGER(CMISSIntg), PARAMETER :: iCaFieldUserNumber=11


  INTEGER(CMISSIntg), PARAMETER :: CaTnCFieldUserNumber=13
  INTEGER(CMISSIntg), PARAMETER :: NumRyRFieldUserNumber=16
  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=17
  INTEGER(CMISSIntg), PARAMETER :: CellMLUserNumber=18
  INTEGER(CMISSIntg), PARAMETER :: CellMLModelsFieldUserNumber=19
  INTEGER(CMISSIntg), PARAMETER :: CellMLStateFieldUserNumber=20
  INTEGER(CMISSIntg), PARAMETER :: CellMLIntermediateFieldUserNumber=21
  INTEGER(CMISSIntg), PARAMETER :: CellMLParametersFieldUserNumber=22
  INTEGER(CMISSIntg), PARAMETER :: GeneratedMeshUserNumber=23

  INTEGER(CMISSIntg), PARAMETER :: FCaEquationsSetUserNumber=24
  INTEGER(CMISSIntg), PARAMETER :: FCaMaterialsFieldUserNumber=25
  INTEGER(CMISSIntg), PARAMETER :: FCaFieldUserNumber=26
  INTEGER(CMISSIntg), PARAMETER :: FCaEquationsSetFieldUserNumber=27
  INTEGER(CMISSIntg), PARAMETER :: iFCaFieldUserNumber=28

  INTEGER(CMISSIntg), PARAMETER :: FEquationsSetUserNumber=29
  INTEGER(CMISSIntg), PARAMETER :: FMaterialsFieldUserNumber=30
  INTEGER(CMISSIntg), PARAMETER :: FFieldUserNumber=31
  INTEGER(CMISSIntg), PARAMETER :: FEquationsSetFieldUserNumber=32
  INTEGER(CMISSIntg), PARAMETER :: iFFieldUserNumber=33

  INTEGER(CMISSIntg), PARAMETER :: CaMEquationsSetUserNumber=34
  INTEGER(CMISSIntg), PARAMETER :: CaMMaterialsFieldUserNumber=35
  INTEGER(CMISSIntg), PARAMETER :: CaMFieldUserNumber=36
  INTEGER(CMISSIntg), PARAMETER :: CaMEquationsSetFieldUserNumber=37
  INTEGER(CMISSIntg), PARAMETER :: iCaMFieldUserNumber=38

  INTEGER(CMISSIntg), PARAMETER :: CaMCaEquationsSetUserNumber=39
  INTEGER(CMISSIntg), PARAMETER :: CaMCaMaterialsFieldUserNumber=40
  INTEGER(CMISSIntg), PARAMETER :: CaMCaFieldUserNumber=41
  INTEGER(CMISSIntg), PARAMETER :: CaMCaEquationsSetFieldUserNumber=42
  INTEGER(CMISSIntg), PARAMETER :: iCaMCaFieldUserNumber=43

  INTEGER(CMISSIntg), PARAMETER :: ATPEquationsSetUserNumber=44
  INTEGER(CMISSIntg), PARAMETER :: ATPMaterialsFieldUserNumber=45
  INTEGER(CMISSIntg), PARAMETER :: ATPFieldUserNumber=46
  INTEGER(CMISSIntg), PARAMETER :: ATPEquationsSetFieldUserNumber=47
  INTEGER(CMISSIntg), PARAMETER :: iATPFieldUserNumber=48

  INTEGER(CMISSIntg), PARAMETER :: ATPCaEquationsSetUserNumber=49
  INTEGER(CMISSIntg), PARAMETER :: ATPCaMaterialsFieldUserNumber=50
  INTEGER(CMISSIntg), PARAMETER :: ATPCaFieldUserNumber=51
  INTEGER(CMISSIntg), PARAMETER :: ATPCaEquationsSetFieldUserNumber=52
  INTEGER(CMISSIntg), PARAMETER :: iATPCaFieldUserNumber=53
  INTEGER(CMISSIntg), PARAMETER :: CaDyadFieldUserNumber=54
  INTEGER(CMISSIntg), PARAMETER :: CaJSRFieldUserNumber=55
  INTEGER(CMISSIntg), PARAMETER :: PopenFieldUserNumber=56
  INTEGER(CMISSIntg), PARAMETER :: RyRDenseFieldUserNumber=57

  INTEGER(CMISSIntg), PARAMETER :: CaSLBufFieldUserNumber=58
  INTEGER(CMISSIntg), PARAMETER :: CaSRBufFieldUserNumber=59
  INTEGER(CMISSIntg), PARAMETER :: iLCCFieldUserNumber=60
  
  


  !CMFE_ variables
  !setting up fields for Ca and CaM equation sets.
  TYPE(CMFE_BasisType) :: Basis
  TYPE(CMFE_CoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(CMFE_DecompositionType) :: Decomposition
  TYPE(CMFE_FieldType) :: GeometricField,CaMaterialsField,FCaMaterialsField,FMaterialsField,CaField,FCaField,FField
  TYPE(CMFE_FieldType) :: CaMMaterialsField,CaMCaMaterialsField,ATPMaterialsField,ATPCaMaterialsField
  TYPE(CMFE_FieldType) :: ATPCaField,CaMField,CaMCaField,RyRDenseField,ATPField,PopenField,CaDyadField,CaJSRField  
  TYPE(CMFE_FieldType) :: CaEquationsSetField,FCaEquationsSetField,FEquationsSetField
  TYPE(CMFE_FieldType) :: CaMEquationsSetField,CaMCaEquationsSetField,ATPEquationsSetField,ATPCaEquationsSetField   
  TYPE(CMFE_FieldsType) :: Fields
  TYPE(CMFE_MeshType) :: Mesh
  TYPE(CMFE_MeshElementsType) :: MeshElements
  TYPE(CMFE_NodesType) :: Nodes
  TYPE(CMFE_RegionType) :: Region,WorldRegion
  TYPE(CMFE_EquationsType) :: CaEquations,FCaEquations,FEquations
  TYPE(CMFE_EquationsType) :: CaMEquations,CaMCaEquations,ATPEquations,ATPCaEquations
  TYPE(CMFE_EquationsSetType) :: CaEquationsSet,FCaEquationsSet,FEquationsSet
  TYPE(CMFE_EquationsSetType) :: CaMEquationsSet,CaMCaEquationsSet,ATPEquationsSet,ATPCaEquationsSet
  TYPE(CMFE_ControlLoopType) :: ControlLoop
  TYPE(CMFE_ProblemType) :: Problem
  TYPE(CMFE_SolverType) :: Solver,LinearSolver
  TYPE(CMFE_SolverEquationsType) :: SolverEquations
  TYPE(CMFE_BoundaryConditionsType) :: BoundaryConditions
  TYPE(CMFE_CellMLType) :: CellML
  TYPE(CMFE_CellMLEquationsType) :: CellMLEquations
  TYPE(CMFE_FieldType) :: CellMLModelsField,CellMLStateField,CellMLIntermediateField,CellMLParametersField
  TYPE(CMFE_FieldType) :: iCaField,CaTnCField,NumRyRField,iFCaField,iFField,iCaMField,iATPField,iCaMCaField,iATPCaField
  TYPE(CMFE_FieldType) :: CaSLBufField,CaSRBufField,iLCCField
  TYPE(CMFE_GeneratedMeshType) :: GeneratedMesh    

  !Program variables
  INTEGER(CMISSIntg) :: NUMBER_OF_ATTRIBUTES,BOUNDARY_MARKER,ELE_ATTRIBUTES
  INTEGER :: st,i,NUMBER_OF_COORDS,NODES_PER_ELE
  INTEGER :: NUMBER_OF_NODES,node,NUMBER_OF_ELEMENTS,element,NUMBER_OF_RYRS,NODE_BD_LABEL
  REAL(CMISSDP),ALLOCATABLE,DIMENSION(:,:) :: NodeCoords
  REAL(CMISSDP), ALLOCATABLE, DIMENSION(:) :: RyRDensity
  INTEGER(CMISSIntg),ALLOCATABLE,DIMENSION(:,:) :: ElemMap
  INTEGER(CMISSIntg),ALLOCATABLE,DIMENSION(:,:) :: NodeNums
  REAL(CMISSDP) :: nodex,nodey,nodez,sphere_constant,RELEASE_RADIUS
  LOGICAL :: EXPORT_FIELD=.FALSE.
  INTEGER(CMISSIntg) :: CELL_TYPE,NumRyRsPerCluster,RELEASE_NODE_RANK
  INTEGER(CMISSIntg) :: ryrModelIndex,GeometricMeshComponent
  INTEGER(CMISSIntg) :: Err,EquationsSetIndex,NODE_NUMBER,RYR_NODE_NUMBER,CONDITION,CellMLIndex,ELEM_NUMBER,NUMBER_RELEASE_NODES
  INTEGER(CMISSIntg),DIMENSION(166) :: CELLBOUNDARYNODES
  INTEGER(CMISSIntg) :: SL_BD_MARKER,MITO_BD_MARKER,MITO_REGION_MARKER,WITH_MITO_ELEMENTS,ELEM_LABEL
  REAL(CMISSDP) :: startT,endT,Tstep,ODE_TIME_STEP,VALUE,init_Ca, init_FCa, init_F, ryr_nodex,ryr_nodey,ryr_nodez, &
    & caDiffx, caDiffy,caDiffz,fcaDiffx, fcaDiffy,fcaDiffz,fDiffx,fDiffy,fDiffz,store_coeff,iCa,init_CaTnC,NodeRyRDensity
  REAL(CMISSDP) :: init_CaM, init_CaMCa, init_ATP, init_ATPCa, &
    & camDiffx, camDiffy,camDiffz,camcaDiffx, camcaDiffy,camcaDiffz,atpDiffx,atpDiffy,atpDiffz, &
    & atpcaDiffx,atpcaDiffy,atpcaDiffz,init_CaDyad,init_CaJSR,init_CaSLBuf,init_CaSRBuf,init_Popen   
  INTEGER(CMISSIntg) :: NonZeroNodes
  CHARACTER(250) :: CELLID,NODEFILE,ELEMFILE,CELLPATH,RyRModel,RYRDENSITYFILE
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber,NodeDomain,ElementDomain
  INTEGER :: MPI_IERROR
#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
    
#ifdef WIN32
  !Initialise QuickWin
  QUICKWIN_WINDOW_CONFIG%TITLE="General Output" !Window title
  QUICKWIN_WINDOW_CONFIG%NUMTEXTROWS=-1 !Max possible number of rows
  QUICKWIN_WINDOW_CONFIG%MODE=QWIN$SCROLLDOWN
  !Set the window parameters
  QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
  !If attempt fails set with system estimated values
  IF(.NOT.QUICKWIN_STATUS) QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
#endif


!_________________________________________________________________________________________________
  !Problem INPUTS. PARAMETERS FROM SOELLER et.al. 2009
  !MESH FILES
  open(unit=9,file='inputs.txt',status='old',action='read',iostat=st)
  IF(st>0)then
    print *,'Error opening inputs file',st
    STOP
  ELSE
    PRINT *,'inputs file opened correctly'
    READ(9,*) !file info
    READ(9,*) !comment on next set of variables - reading in file info
    READ(9,*) RYR_NODE_NUMBER
    READ(9,*)
    READ(9,*) NumRyRsPerCluster,NodeRyRDensity
    READ(9,*)
    READ(9,*) RyRModel
    READ(9,*)
    READ(9,*) init_Ca,caDiffx,caDiffy,caDiffz,init_Popen
    READ(9,*)
    READ(9,*) init_F,fDiffx,fDiffy,fDiffz
    READ(9,*)
    READ(9,*) init_FCa,fcaDiffx,fcaDiffy,fcaDiffz
    READ(9,*)
    READ(9,*) init_CaM,camDiffx,camDiffy,camDiffz
    READ(9,*)
    READ(9,*) init_CaMCa,camcaDiffx,camcaDiffy,camcaDiffz
    READ(9,*)
    READ(9,*) init_ATP,atpDiffx,atpDiffy,atpDiffz
    READ(9,*)
    READ(9,*) init_ATPCa,atpcaDiffx,atpcaDiffy,atpcaDiffz
    READ(9,*)
    READ(9,*) init_CaTnC,init_CaSRBuf,init_CaSLBuf,init_CaDyad,init_CaJSR
    READ(9,*) 
    READ(9,*) startT,endT,Tstep,ODE_TIME_STEP
  ENDIF
  CLOSE(9)
  !Write the params out to screen for double checking
  WRITE(*,*) 'RyR Cluster Node:', RYR_NODE_NUMBER
  WRITE(*,*) 'Number of RyRs per Cluster:', NumRyRsPerCluster
  WRITE(*,*) 'RyR Clusters per Unit Volume At Node:', NodeRyRDensity
  WRITE(*,*) 'CellML Model File:', RyRModel
  !cell initial conditions
  WRITE(*,*) 'Initial [Ca]i = ',init_Ca !dependent field (cytosolic calcium) set to 0.1 microM at rest.
  WRITE(*,*) 'Ca Diff Coeff in x = ',caDiffx
  WRITE(*,*) 'Ca Diff Coeff in y = ',caDiffy
  WRITE(*,*) 'Ca Diff Coeff in z = ',caDiffz

  WRITE(*,*) 'Initial [F]i = ',init_F !dependent field (cytosolic calcium) set to 0.1 microM at rest.
  WRITE(*,*) 'F Diff Coeff in x = ',fDiffx
  WRITE(*,*) 'F Diff Coeff in y = ',fDiffy
  WRITE(*,*) 'F Diff Coeff in z = ',fDiffz

  WRITE(*,*) 'Initial [FCa]i = ',init_FCa !dependent field (cytosolic calcium) set to 0.1 microM at rest.
  WRITE(*,*) 'FCa Diff Coeff in x = ',fcaDiffx
  WRITE(*,*) 'FCa Diff Coeff in y = ',fcaDiffy
  WRITE(*,*) 'FCa Diff Coeff in z = ',fcaDiffz

  WRITE(*,*) 'Initial Equil. [CaTnC] - uniform over myofibril region = ',init_CaTnC
  store_coeff = 1.0_CMISSDP

  WRITE(*,*) 'Tstart=',startT
  WRITE(*,*) 'Tend=',endT
  WRITE(*,*) 'Tstep=',Tstep
  WRITE(*,*) 'ODE_Tstep=',ODE_TIME_STEP
  
  EXPORT_FIELD=.FALSE.
!_________________________________________________________________________________________________
  !Intialise OpenCMFE_
  CALL CMFE_Initialise(WorldCoordinateSystem,WorldRegion,Err)
  CALL CMFE_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,Err)
  !get computational nodes for parallel processing
  CALL CMFE_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL CMFE_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)

   !Start the creation of a new RC coordinate system
  CALL CMFE_CoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL CMFE_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  CALL CMFE_CoordinateSystem_DimensionSet(CoordinateSystemUserNumber,3,Err)
  !The coordinate system is 3D by default;set it to be 3D in above command.
  !Finish the creation of the coordinate system
  CALL CMFE_CoordinateSystem_CreateFinish(CoordinateSystem,Err)

  !Start the creation of the region
  CALL CMFE_Region_Initialise(Region,Err)
  CALL CMFE_Region_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  !Set the regions coordinate system to the 3D RC coordinate system that we have created
  CALL CMFE_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
  CALL CMFE_Region_LabelSet(Region,"Cube",Err)
  !Finish the creation of the region
  CALL CMFE_Region_CreateFinish(Region,Err)

!_________________________________________________________________________________________________
  !Start the creation of a trilinear-simplex basis
  CALL CMFE_Basis_Initialise(Basis,Err)
  CALL CMFE_Basis_CreateStart(BasisUserNumber,Basis,Err)
  !Set the basis to be a trilinear simplex  basis
  CALL CMFE_Basis_TypeSet(Basis,CMFE_BASIS_SIMPLEX_TYPE,Err)
  CALL CMFE_Basis_NumberOfXiSet(Basis,3,Err)
  !set interpolation to be linear
  CALL CMFE_Basis_InterpolationXiSet(Basis,(/CMFE_Basis_Linear_Simplex_Interpolation, &
   &   CMFE_Basis_Linear_Simplex_Interpolation, CMFE_Basis_Linear_Simplex_Interpolation/),Err)
  !Finish the creation of the basis
  CALL CMFE_Basis_CreateFinish(Basis,Err)

!_________________________________________________________________________________________________  
  !Time to create a mesh - wohoo!
  !Read in nodes (set up RyRDensity array with column 
  !of zeros for later updating).
  NODEFILE="cuboid.1.node"
  ELEMFILE="cuboid.1.ele"
  open(unit=10,file=NODEFILE,status='old',action='read',iostat=st)
  IF(st>0)then
    print *,'Error opening node file',st
    STOP
  ELSE
    PRINT *,'Node file opened correctly'
    READ(10,*) NUMBER_OF_NODES, NUMBER_OF_COORDS, NUMBER_OF_ATTRIBUTES, BOUNDARY_MARKER
    ALLOCATE(NodeNums(NUMBER_OF_NODES,2))
    ALLOCATE(NodeCoords(NUMBER_OF_NODES,NUMBER_OF_COORDS))
    DO i = 1,NUMBER_OF_NODES
      READ(10,*) NodeNums(i,1),NodeCoords(i,1),NodeCoords(i,2),NodeCoords(i,3),NodeNums(i,2)
    ENDDO
  ENDIF
  CLOSE(10)
  !Read in elements
  OPEN(unit=11,file=ELEMFILE,status='old',action='read',iostat=st)
  IF(st>0)THEN
    PRINT *,'Error opening element file',st
    STOP
  ELSE
    PRINT *,'Element file opened successfully'
    READ(11,*) NUMBER_OF_ELEMENTS,NODES_PER_ELE,ELE_ATTRIBUTES
    ALLOCATE(ElemMap(NUMBER_OF_ELEMENTS,5))
    DO i = 1,NUMBER_OF_ELEMENTS
      READ(11,*) ElemMap(i,1),ElemMap(i,2),ElemMap(i,3),ElemMap(i,4),ElemMap(i,5)
    ENDDO
  ENDIF 
  CLOSE(11)

  CALL CMFE_Nodes_Initialise(Nodes,Err)
  CALL CMFE_Nodes_CreateStart(Region,NUMBER_OF_NODES,Nodes,Err)
  CALL CMFE_Nodes_CreateFinish(Nodes,Err)
  
  CALL CMFE_Mesh_Initialise(Mesh,Err)
  CALL CMFE_Mesh_CreateStart(MeshUserNumber,Region,NUMBER_OF_COORDS,Mesh,Err)
  CALL CMFE_Mesh_NumberOfElementsSet(Mesh,NUMBER_OF_ELEMENTS,Err)
  CALL CMFE_Mesh_NumberOfComponentsSet(Mesh,1,Err)
  
  CALL CMFE_MeshElements_Initialise(MeshElements,Err)
  CALL CMFE_MeshElements_CreateStart(Mesh,1,Basis,MeshElements,Err)
  DO i = 1,NUMBER_OF_ELEMENTS
    element = ElemMap(i,1)
    CALL CMFE_MeshElements_NodesSet(MeshElements,element,(/ElemMap(i,2),ElemMap(i,3), &
     &   ElemMap(i,4),ElemMap(i,5)/),Err)
  ENDDO
  CALL CMFE_MeshElements_CreateFinish(MeshElements,Err)
  CALL CMFE_Mesh_CreateFinish(Mesh,Err)

  !Create a decomposition
  CALL CMFE_Decomposition_Initialise(Decomposition,Err)
  CALL CMFE_Decomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL CMFE_Decomposition_TypeSet(Decomposition,CMFE_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL CMFE_Decomposition_NumberOfDomainsSet(Decomposition,NumberOfComputationalNodes,Err)
  !Finish the decomposition
  CALL CMFE_Decomposition_CreateFinish(Decomposition,Err)

  !Start to create a default (geometric) field on the region
  CALL CMFE_Field_Initialise(GeometricField,Err)
  CALL CMFE_Field_CreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  !Set the decomposition to use
  CALL CMFE_Field_MeshDecompositionSet(GeometricField,Decomposition,Err)
  !Set the domain to be used by the field components. We have 3 field components in 1 mesh component
  CALL CMFE_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL CMFE_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,2,1,Err)
  CALL CMFE_Field_ComponentMeshComponentSet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,3,1,Err)
  !Finish creating the field
  CALL CMFE_Field_CreateFinish(GeometricField,Err)

  !Set the geometric field values

  DO i = 1,NUMBER_OF_NODES
    node = NodeNums(i,1)
    CALL CMFE_Decomposition_NodeDomainGet(Decomposition,node,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      nodex = NodeCoords(i,1)
      nodey = NodeCoords(i,2)
      nodez = NodeCoords(i,3)
      CALL CMFE_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
       &   1,node,1,nodex,Err)
      CALL CMFE_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
       &   1,node,2,nodey,Err)
      CALL CMFE_Field_ParameterSetUpdateNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
       &   1,node,3,nodez,Err)
     ENDIF
    ENDDO
  CALL CMFE_Field_ParameterSetUpdateStart(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL CMFE_Field_ParameterSetUpdateFinish(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)

  IF(EXPORT_FIELD) THEN
    CALL CMFE_Fields_Initialise(Fields,Err)
    CALL CMFE_Fields_Create(Region,Fields,Err)
    CALL CMFE_Fields_NodesExport(Fields,"Cube_Geom","FORTRAN",Err)
    CALL CMFE_Fields_ElementsExport(Fields,"Cube_Geom","FORTRAN",Err)
    CALL CMFE_Fields_Finalise(Fields,Err)
  ENDIF 

!______________________________________________________________________________________________________________
  WRITE(*,*) 'Create the cellml reaction with split reaction diffusion equations_set - 1 for each species' 
  !Ca equations
  CALL CMFE_EquationsSet_Initialise(CaEquationsSet,Err)
  CALL CMFE_Field_Initialise(CaEquationsSetField,Err)
  CALL CMFE_EquationsSet_CreateStart(CaEquationsSetUserNumber,Region, & 
    & GeometricField,[CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMFE_EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE, &
    & CMFE_EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE], &
    & CaEquationsSetFieldUserNumber,CaEquationsSetField,CaEquationsSet,Err)
  !Set the equations set to be a standard Diffusion no source problem
  !Finish creating the equations set
  CALL CMFE_EquationsSet_CreateFinish(CaEquationsSet,Err)


  !Create the equations set dependent field variables for Ca
  CALL CMFE_Field_Initialise(CaField,Err)
  CALL CMFE_EquationsSet_DependentCreateStart(CaEquationsSet,CaFieldUserNumber,CaField,Err)
  CALL CMFE_Field_VariableLabelSet(CaField,CMFE_FIELD_U_VARIABLE_TYPE,"Ca Field",Err)
  !Finish the equations set dependent field variables
  CALL CMFE_EquationsSet_DependentCreateFinish(CaEquationsSet,Err)
  CALL CMFE_Field_ComponentValuesInitialise(CaField,CMFE_FIELD_U_VARIABLE_TYPE, &
    & CMFE_FIELD_VALUES_SET_TYPE,1,init_Ca,Err)


  !Create the equations set material field variables - Ca
  !by default 2 comps for reac diff i.e. diff coeff in 1 direction set constant spatially = 1, and storage coeff set to 1
  CALL CMFE_Field_Initialise(CaMaterialsField,Err)
  CALL CMFE_EquationsSet_MaterialsCreateStart(CaEquationsSet,CaMaterialsFieldUserNumber,CaMaterialsField,Err)
  CALL CMFE_Field_VariableLabelSet(CaMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,"Ca Materials Field",Err)
  !Finish the equations set materials field variables
  CALL CMFE_EquationsSet_MaterialsCreateFinish(CaEquationsSet,Err)
  CALL CMFE_Field_ComponentValuesInitialise(CaMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
   & 1,caDiffx,Err) !ca diff coeff in x
  CALL CMFE_Field_ComponentValuesInitialise(CaMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
   & 2,caDiffy,Err) !ca diff coeff in y
  CALL CMFE_Field_ComponentValuesInitialise(CaMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
   & 3,caDiffz,Err) !ca diff coeff in z
  CALL CMFE_Field_ComponentValuesInitialise(CaMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
   & 4,store_coeff,Err) ! storage coefficient

  CALL CMFE_Field_ParameterSetUpdateStart(CaMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL CMFE_Field_ParameterSetUpdateFinish(CaMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)

 
  !Set up source field for reaction diffusion equation set. Note that for the split problem subtype, the source field is not used at all.
  !iCaField
  !Might use the field for CellML input of elementary RyR calcium release
  CALL CMFE_Field_Initialise(PopenField,Err)
  CALL CMFE_EquationsSet_SourceCreateStart(CaEquationsSet,PopenFieldUserNumber,PopenField,Err)
  CALL CMFE_Field_VariableLabelSet(PopenField,CMFE_FIELD_U_VARIABLE_TYPE,"Popen Field",Err)
  !Finish the equations set source field variables
  CALL CMFE_EquationsSet_SourceCreateFinish(CaEquationsSet,Err)
  !Initialising the iCaField to iCa everywhere. Modifying for RyRs in a later loop.
  CALL CMFE_Field_ComponentValuesInitialise(PopenField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,init_Popen,Err)

  !CaSLBuf
  CALL CMFE_Field_Initialise(CaSLBufField,Err)
  CALL CMFE_Field_CreateStart(CaSLBufFieldUserNumber,Region,CaSLBufField,Err)
  CALL CMFE_Field_TypeSet(CaSLBufField,CMFE_FIELD_GENERAL_TYPE,Err)
  CALL CMFE_Field_MeshDecompositionSet(CaSLBufField,Decomposition,Err)
  CALL CMFE_Field_GeometricFieldSet(CaSLBufField,GeometricField,Err)
  CALL CMFE_Field_NumberOfVariablesSet(CaSLBufField,1,Err)
  CALL CMFE_Field_VariableTypesSet(CaSLBufField,[CMFE_FIELD_U_VARIABLE_TYPE],Err)
  CALL CMFE_Field_DataTypeSet(CaSLBufField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_DP_TYPE,Err)
  CALL CMFE_Field_DimensionSet(CaSLBufField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_SCALAR_DIMENSION_TYPE,Err)
  CALL CMFE_Field_NumberOfComponentsSet(CaSLBufField,CMFE_FIELD_U_VARIABLE_TYPE,1,Err)
  CALL CMFE_Field_VariableLabelSet(CaSLBufField,CMFE_FIELD_U_VARIABLE_TYPE,"CaSLBuf Field",Err)
  CALL CMFE_Field_ComponentMeshComponentGet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE, & 
    & 1,GeometricMeshComponent,ERR)
  !Default to the geometric interpolation setup
  CALL CMFE_Field_ComponentMeshComponentSet(CaSLBufField,CMFE_FIELD_U_VARIABLE_TYPE,1, &
    & GeometricMeshComponent,ERR)            
  !Specify the interpolation to be same as geometric interpolation
  CALL CMFE_Field_ComponentInterpolationSet(CaSLBufField,CMFE_FIELD_U_VARIABLE_TYPE,1, &
    & CMFE_FIELD_NODE_BASED_INTERPOLATION,ERR)
  CALL CMFE_Field_CreateFinish(CaSLBufField,Err)
  !Initialise CaSLBuf concentration to equilibrium value
  !Set the values to be nodally varying - mito nodes with different concentrations than myo regions
  CALL CMFE_Field_ComponentValuesInitialise(CaSLBufField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,0.0_CMISSDP,Err)
  CALL CMFE_Field_ParameterSetUpdateStart(CaSLBufField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL CMFE_Field_ParameterSetUpdateFinish(CaSLBufField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)


  !iLCCField
  !Set up number of RyRs in each cluster as a spatial distribution field
  !set up intensity field
  CALL CMFE_Field_Initialise(iLCCField,Err)
  CALL CMFE_Field_CreateStart(iLCCFieldUserNumber,Region,iLCCField,Err)
  CALL CMFE_Field_TypeSet(iLCCField,CMFE_FIELD_GENERAL_TYPE,Err)
  CALL CMFE_Field_MeshDecompositionSet(iLCCField,Decomposition,Err)
  CALL CMFE_Field_GeometricFieldSet(iLCCField,GeometricField,Err)
  CALL CMFE_Field_NumberOfVariablesSet(iLCCField,1,Err)
  CALL CMFE_Field_VariableTypesSet(iLCCField,[CMFE_FIELD_U_VARIABLE_TYPE],Err)
  CALL CMFE_Field_DataTypeSet(iLCCField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_DP_TYPE,Err)
  CALL CMFE_Field_DimensionSet(iLCCField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_SCALAR_DIMENSION_TYPE,Err)
  CALL CMFE_Field_NumberOfComponentsSet(iLCCField,CMFE_FIELD_U_VARIABLE_TYPE,1,Err)
  CALL CMFE_Field_VariableLabelSet(iLCCField,CMFE_FIELD_U_VARIABLE_TYPE,"iLCC Field",Err)
  CALL CMFE_Field_ComponentMeshComponentGet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE, & 
    & 1,GeometricMeshComponent,ERR)
  !Default to the geometric interpolation setup
  CALL CMFE_Field_ComponentMeshComponentSet(iLCCField,CMFE_FIELD_U_VARIABLE_TYPE,1, &
    & GeometricMeshComponent,ERR)            
  !Specify the interpolation to be same as geometric interpolation
  CALL CMFE_Field_ComponentInterpolationSet(iLCCField,CMFE_FIELD_U_VARIABLE_TYPE,1, &
    & CMFE_FIELD_NODE_BASED_INTERPOLATION,ERR)
  CALL CMFE_Field_CreateFinish(iLCCField,Err)
  !Initialise Num RyR Field
  CALL CMFE_Field_ComponentValuesInitialise(iLCCField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,0.0_CMISSDP,Err)


  !NumRyRField
  !Set up number of RyRs in each cluster as a spatial distribution field
  !set up intensity field
  CALL CMFE_Field_Initialise(NumRyRField,Err)
  CALL CMFE_Field_CreateStart(NumRyRFieldUserNumber,Region,NumRyRField,Err)
  CALL CMFE_Field_TypeSet(NumRyRField,CMFE_FIELD_GENERAL_TYPE,Err)
  CALL CMFE_Field_MeshDecompositionSet(NumRyRField,Decomposition,Err)
  CALL CMFE_Field_GeometricFieldSet(NumRyRField,GeometricField,Err)
  CALL CMFE_Field_NumberOfVariablesSet(NumRyRField,1,Err)
  CALL CMFE_Field_VariableTypesSet(NumRyRField,[CMFE_FIELD_U_VARIABLE_TYPE],Err)
  CALL CMFE_Field_DataTypeSet(NumRyRField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_DP_TYPE,Err)
  CALL CMFE_Field_DimensionSet(NumRyRField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_SCALAR_DIMENSION_TYPE,Err)
  CALL CMFE_Field_NumberOfComponentsSet(NumRyRField,CMFE_FIELD_U_VARIABLE_TYPE,1,Err)
  CALL CMFE_Field_VariableLabelSet(NumRyRField,CMFE_FIELD_U_VARIABLE_TYPE,"Num RyR Field",Err)
  CALL CMFE_Field_ComponentMeshComponentGet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE, & 
    & 1,GeometricMeshComponent,ERR)
  !Default to the geometric interpolation setup
  CALL CMFE_Field_ComponentMeshComponentSet(NumRyRField,CMFE_FIELD_U_VARIABLE_TYPE,1, &
    & GeometricMeshComponent,ERR)            
  !Specify the interpolation to be same as geometric interpolation
  CALL CMFE_Field_ComponentInterpolationSet(NumRyRField,CMFE_FIELD_U_VARIABLE_TYPE,1, &
    & CMFE_FIELD_NODE_BASED_INTERPOLATION,ERR)
  CALL CMFE_Field_CreateFinish(NumRyRField,Err)
  !Initialise Num RyR Field
  CALL CMFE_Field_ComponentValuesInitialise(NumRyRField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,0.0_CMISSDP,Err)

  !Define a single RyR release node, multiply the intensity by number of ryrs per cluster
  CALL CMFE_Decomposition_NodeDomainGet(Decomposition,RYR_NODE_NUMBER,1,NodeDomain,Err)
  RELEASE_NODE_RANK = NodeDomain

  IF(NodeDomain.EQ.ComputationalNodeNumber) THEN
    RELEASE_NODE_RANK = ComputationalNodeNumber
    CALL CMFE_Field_ParameterSetUpdateNode(NumRyRField,CMFE_FIELD_U_VARIABLE_TYPE, &
      & CMFE_FIELD_VALUES_SET_TYPE,1,1,RYR_NODE_NUMBER,1,(NumRyRsPerCluster*NodeRyRDensity),Err)

    !Setting SLBuffer concentration only at nodes that make the dyad.
    CALL CMFE_Field_ParameterSetUpdateNode(CaSLBufField,CMFE_FIELD_U_VARIABLE_TYPE, &
      & CMFE_FIELD_VALUES_SET_TYPE,1,1,RYR_NODE_NUMBER,1,(init_CaSLBuf),Err)

    CALL CMFE_Field_ParameterSetUpdateNode(iLCCField,CMFE_FIELD_U_VARIABLE_TYPE, &
      & CMFE_FIELD_VALUES_SET_TYPE,1,1,RYR_NODE_NUMBER,1,0.5e-16_CMISSDP,Err)

    !noting the coordinates of the ryr release node.
    CALL CMFE_Field_ParameterSetGetNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
      &   1,RYR_NODE_NUMBER,1,ryr_nodex,Err)
    CALL CMFE_Field_ParameterSetGetNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
      &   1,RYR_NODE_NUMBER,2,ryr_nodey,Err)
    CALL CMFE_Field_ParameterSetGetNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
      &   1,RYR_NODE_NUMBER,3,ryr_nodez,Err)

  ENDIF
  !CALL MPI_BARRIER(MPI_COMM_WORLD,MPI_IERROR)

  !CALL MPI_BCAST(ryr_nodex,1,MPI_DOUBLE,RELEASE_NODE_RANK,MPI_COMM_WORLD,MPI_IERROR)
  !CALL MPI_BCAST(ryr_nodey,1,MPI_DOUBLE,RELEASE_NODE_RANK,MPI_COMM_WORLD,MPI_IERROR)
  !CALL MPI_BCAST(ryr_nodez,1,MPI_DOUBLE,RELEASE_NODE_RANK,MPI_COMM_WORLD,MPI_IERROR)

  CALL CMFE_Field_ParameterSetUpdateStart(NumRyRField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL CMFE_Field_ParameterSetUpdateFinish(NumRyRField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL CMFE_Field_ParameterSetUpdateStart(CaSLBufField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL CMFE_Field_ParameterSetUpdateFinish(CaSLBufField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  
  !Now assign ryr density to nodes that fall within a 100 nm (0.1 micron) radius of the center of the RyR_Node_Number
  RELEASE_RADIUS = 0.1_CMISSDP
  NUMBER_RELEASE_NODES=1
  DO node=1,NUMBER_OF_NODES
    NODE_NUMBER=NodeNums(node,1)
    CALL CMFE_Decomposition_NodeDomainGet(Decomposition,NODE_NUMBER,1,NodeDomain,Err)
    IF(NodeDomain.EQ.ComputationalNodeNumber) THEN
      CALL CMFE_Field_ParameterSetGetNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
        &   1,NODE_NUMBER,1,nodex,Err)
      CALL CMFE_Field_ParameterSetGetNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
        &   1,NODE_NUMBER,2,nodey,Err)
      CALL CMFE_Field_ParameterSetGetNode(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1, &
        &   1,NODE_NUMBER,3,nodez,Err)
      sphere_constant = SQRT((ryr_nodex-nodex)**2+(ryr_nodey-nodey)**2+(ryr_nodez-nodez)**2)
      IF(sphere_constant.LE.RELEASE_RADIUS) THEN
        WRITE(*,*) 'Node Number',NODE_NUMBER
        WRITE(*,*) 'COORDS:', nodex,nodey,nodez
        WRITE(*,*) 'REL CENTRE:',ryr_nodex,ryr_nodey,ryr_nodez
        WRITE(*,*) 'DISTANCE:',sphere_constant 
        CALL CMFE_Field_ParameterSetUpdateNode(NumRyRField,CMFE_FIELD_U_VARIABLE_TYPE, &
          & CMFE_FIELD_VALUES_SET_TYPE,1,1,NODE_NUMBER,1,(NumRyRsPerCluster*NodeRyRDensity),Err)
            
        !Setting SLBuffer concentration only at nodes that make the dyad.
        CALL CMFE_Field_ParameterSetUpdateNode(CaSLBufField,CMFE_FIELD_U_VARIABLE_TYPE, &
          & CMFE_FIELD_VALUES_SET_TYPE,1,1,NODE_NUMBER,1,(init_CaSLBuf),Err)

        CALL CMFE_Field_ParameterSetUpdateNode(iLCCField,CMFE_FIELD_U_VARIABLE_TYPE, &
          & CMFE_FIELD_VALUES_SET_TYPE,1,1,NODE_NUMBER,1,0.5e-16_CMISSDP,Err)

        NUMBER_RELEASE_NODES = NUMBER_RELEASE_NODES+1
      ENDIF
    ENDIF
  ENDDO
  CALL CMFE_Field_ParameterSetUpdateStart(NumRyRField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL CMFE_Field_ParameterSetUpdateFinish(NumRyRField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL CMFE_Field_ParameterSetUpdateStart(iLCCField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL CMFE_Field_ParameterSetUpdateFinish(iLCCField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL CMFE_Field_ParameterSetUpdateStart(CaSLBufField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL CMFE_Field_ParameterSetUpdateFinish(CaSLBufField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  
  WRITE(*,*) "The final number of nodes which will be releasing Ca2+ =", NUMBER_RELEASE_NODES
  !Set up the fields for the other buffers which will store concentrations of the Ca-Buffer complex
  !F equations
  CALL CMFE_EquationsSet_Initialise(FEquationsSet,Err)
  CALL CMFE_Field_Initialise(FEquationsSetField,Err)
  CALL CMFE_EquationsSet_CreateStart(FEquationsSetUserNumber,Region, & 
    & GeometricField,[CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMFE_EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE, &
    & CMFE_EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE], &
    & FEquationsSetFieldUserNumber,FEquationsSetField,FEquationsSet,Err)
  !Set the equations set to be a standard Diffusion no source problem
  !Finish creating the equations set
  CALL CMFE_EquationsSet_CreateFinish(FEquationsSet,Err)


  !Create the equations set dependent field variables for F
  CALL CMFE_Field_Initialise(FField,Err)
  CALL CMFE_EquationsSet_DependentCreateStart(FEquationsSet,FFieldUserNumber,FField,Err)
  CALL CMFE_Field_VariableLabelSet(FField,CMFE_FIELD_U_VARIABLE_TYPE,"F Field",Err)
  !Finish the equations set dependent field variables
  CALL CMFE_EquationsSet_DependentCreateFinish(FEquationsSet,Err)
  CALL CMFE_Field_ComponentValuesInitialise(FField,CMFE_FIELD_U_VARIABLE_TYPE, &
    & CMFE_FIELD_VALUES_SET_TYPE,1,init_F,Err)


  !Create the equations set material field variables - F
  !by default 2 comps for reac diff i.e. diff coeff in 1 direction set constant spatially = 1, and storage coeff set to 1
  CALL CMFE_Field_Initialise(FMaterialsField,Err)
  CALL CMFE_EquationsSet_MaterialsCreateStart(FEquationsSet,FMaterialsFieldUserNumber,FMaterialsField,Err)
  CALL CMFE_Field_VariableLabelSet(FMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,"F Materials Field",Err)
  !Finish the equations set materials field variables
  CALL CMFE_EquationsSet_MaterialsCreateFinish(FEquationsSet,Err)
  CALL CMFE_Field_ComponentValuesInitialise(FMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
   & 1,fDiffx,Err) !f diff coeff in x
  CALL CMFE_Field_ComponentValuesInitialise(FMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
   & 2,fDiffy,Err) !f diff coeff in y
  CALL CMFE_Field_ComponentValuesInitialise(FMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
   & 3,fDiffz,Err) !f diff coeff in z
  CALL CMFE_Field_ComponentValuesInitialise(FMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
   & 4,store_coeff,Err) ! storage coefficient

  CALL CMFE_Field_ParameterSetUpdateStart(FMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL CMFE_Field_ParameterSetUpdateFinish(FMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)

  !Set up source field for reaction diffusion equation set. Note that for the split problem subtype, the source field is not used at all.
  !iFField
  CALL CMFE_Field_Initialise(iFField,Err)
  CALL CMFE_EquationsSet_SourceCreateStart(FEquationsSet,iFFieldUserNumber,iFField,Err)
  CALL CMFE_Field_VariableLabelSet(iFField,CMFE_FIELD_U_VARIABLE_TYPE,"iF Field",Err)
  !Finish the equations set source field variables
  CALL CMFE_EquationsSet_SourceCreateFinish(FEquationsSet,Err)
  !Initialising the iField to zero everywhere. Might modify for RyRs in a later loop.
  CALL CMFE_Field_ComponentValuesInitialise(iFField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,0.0_CMISSDP,Err)

  !FCa equations
  CALL CMFE_EquationsSet_Initialise(FCaEquationsSet,Err)
  CALL CMFE_Field_Initialise(FCaEquationsSetField,Err)
  CALL CMFE_EquationsSet_CreateStart(FCaEquationsSetUserNumber,Region, & 
    & GeometricField,[CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMFE_EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE, &
    & CMFE_EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE], &
    & FCaEquationsSetFieldUserNumber,FCaEquationsSetField,FCaEquationsSet,Err)
  !Set the equations set to be a standard Diffusion no source problem
  !Finish creating the equations set
  CALL CMFE_EquationsSet_CreateFinish(FCaEquationsSet,Err)


  !Create the equations set dependent field variables for FCa
  CALL CMFE_Field_Initialise(FCaField,Err)
  CALL CMFE_EquationsSet_DependentCreateStart(FCaEquationsSet,FCaFieldUserNumber,FCaField,Err)
  CALL CMFE_Field_VariableLabelSet(FCaField,CMFE_FIELD_U_VARIABLE_TYPE,"FCa Field",Err)
  !Finish the equations set dependent field variables
  CALL CMFE_EquationsSet_DependentCreateFinish(FCaEquationsSet,Err)
  CALL CMFE_Field_ComponentValuesInitialise(FCaField,CMFE_FIELD_U_VARIABLE_TYPE, &
    & CMFE_FIELD_VALUES_SET_TYPE,1,init_FCa,Err)


  !Create the equations set material field variables - FCa
  !by default 2 comps for reac diff i.e. diff coeff in 1 direction set constant spatially = 1, and storage coeff set to 1
  CALL CMFE_Field_Initialise(FCaMaterialsField,Err)
  CALL CMFE_EquationsSet_MaterialsCreateStart(FCaEquationsSet,FCaMaterialsFieldUserNumber,FCaMaterialsField,Err)
  CALL CMFE_Field_VariableLabelSet(FCaMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,"FCa Materials Field",Err)
  !Finish the equations set materials field variables
  CALL CMFE_EquationsSet_MaterialsCreateFinish(FCaEquationsSet,Err)
  CALL CMFE_Field_ComponentValuesInitialise(FCaMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
   & 1,fcaDiffx,Err) !fca diff coeff in x
  CALL CMFE_Field_ComponentValuesInitialise(FCaMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
   & 2,fcaDiffy,Err) !fca diff coeff in y
  CALL CMFE_Field_ComponentValuesInitialise(FCaMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
   & 3,fcaDiffz,Err) !fca diff coeff in z
  CALL CMFE_Field_ComponentValuesInitialise(FCaMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
   & 4,store_coeff,Err) ! storage coefficient

  CALL CMFE_Field_ParameterSetUpdateStart(FCaMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL CMFE_Field_ParameterSetUpdateFinish(FCaMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)

  !Set up source field for reaction diffusion equation set. Note that for the split problem subtype, the source field is not used at all.
  !iFCaField
  CALL CMFE_Field_Initialise(iFCaField,Err)
  CALL CMFE_EquationsSet_SourceCreateStart(FCaEquationsSet,iFCaFieldUserNumber,iFCaField,Err)
  CALL CMFE_Field_VariableLabelSet(iFCaField,CMFE_FIELD_U_VARIABLE_TYPE,"iFCa Field",Err)
  !Finish the equations set source field variables
  CALL CMFE_EquationsSet_SourceCreateFinish(FCaEquationsSet,Err)
  !Initialising the iCaField to zero everywhere. Modifying for RyRs in a later loop.
  CALL CMFE_Field_ComponentValuesInitialise(iFCaField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,0.0_CMISSDP,Err)

  !CaM equations
  CALL CMFE_EquationsSet_Initialise(CaMEquationsSet,Err)
  CALL CMFE_Field_Initialise(CaMEquationsSetField,Err)
  CALL CMFE_EquationsSet_CreateStart(CaMEquationsSetUserNumber,Region, & 
    & GeometricField,[CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMFE_EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE, &
    & CMFE_EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE], &
    & CaMEquationsSetFieldUserNumber,CaMEquationsSetField,CaMEquationsSet,Err)
  !Set the equations set to be a standard Diffusion no source problem
  !Finish creating the equations set
  CALL CMFE_EquationsSet_CreateFinish(CaMEquationsSet,Err)


  !Create the equations set dependent field variables for CaM
  CALL CMFE_Field_Initialise(CaMField,Err)
  CALL CMFE_EquationsSet_DependentCreateStart(CaMEquationsSet,CaMFieldUserNumber,CaMField,Err)
  CALL CMFE_Field_VariableLabelSet(CaMField,CMFE_FIELD_U_VARIABLE_TYPE,"CaM Field",Err)
  !Finish the equations set dependent field variables
  CALL CMFE_EquationsSet_DependentCreateFinish(CaMEquationsSet,Err)
  CALL CMFE_Field_ComponentValuesInitialise(CaMField,CMFE_FIELD_U_VARIABLE_TYPE, &
    & CMFE_FIELD_VALUES_SET_TYPE,1,init_CaM,Err)


  !Create the equations set material field variables - CaM
  !by default 2 comps for reac diff i.e. diff coeff in 1 direction set constant spatially = 1, and storage coeff set to 1
  CALL CMFE_Field_Initialise(CaMMaterialsField,Err)
  CALL CMFE_EquationsSet_MaterialsCreateStart(CaMEquationsSet,CaMMaterialsFieldUserNumber,CaMMaterialsField,Err)
  CALL CMFE_Field_VariableLabelSet(CaMMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,"CaM Materials Field",Err)
  !Finish the equations set materials field variables
  CALL CMFE_EquationsSet_MaterialsCreateFinish(CaMEquationsSet,Err)
  CALL CMFE_Field_ComponentValuesInitialise(CaMMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
   & 1,camDiffx,Err) !CaM diff coeff in x
  CALL CMFE_Field_ComponentValuesInitialise(CaMMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
   & 2,camDiffy,Err) !CaM diff coeff in y
  CALL CMFE_Field_ComponentValuesInitialise(CaMMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
   & 3,camDiffz,Err) !CaM diff coeff in z
  CALL CMFE_Field_ComponentValuesInitialise(CaMMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
   & 4,store_coeff,Err) ! storage coefficient

  CALL CMFE_Field_ParameterSetUpdateStart(CaMMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL CMFE_Field_ParameterSetUpdateFinish(CaMMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)

  !Set up source field for reaction diffusion equation set. Note that for the split problem subtype, the source field is not used at all.
  !iCaMField
  CALL CMFE_Field_Initialise(iCaMField,Err)
  CALL CMFE_EquationsSet_SourceCreateStart(CaMEquationsSet,iCaMFieldUserNumber,iCaMField,Err)
  CALL CMFE_Field_VariableLabelSet(iCaMField,CMFE_FIELD_U_VARIABLE_TYPE,"iCaM Field",Err)
  !Finish the equations set source field variables
  CALL CMFE_EquationsSet_SourceCreateFinish(CaMEquationsSet,Err)
  !Initialising the iCaField to zero everywhere. Modifying for RyRs in a later loop.
  CALL CMFE_Field_ComponentValuesInitialise(iCaMField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,0.0_CMISSDP,Err)
    
  !CaMCa equations
  CALL CMFE_EquationsSet_Initialise(CaMCaEquationsSet,Err)
  CALL CMFE_Field_Initialise(CaMCaEquationsSetField,Err)
  CALL CMFE_EquationsSet_CreateStart(CaMCaEquationsSetUserNumber,Region, & 
    & GeometricField,[CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMFE_EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE, &
    & CMFE_EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE], &
    & CaMCaEquationsSetFieldUserNumber,CaMCaEquationsSetField,CaMCaEquationsSet,Err)
  !Set the equations set to be a standard Diffusion no source problem
  !Finish creating the equations set
  CALL CMFE_EquationsSet_CreateFinish(CaMCaEquationsSet,Err)


  !Create the equations set dependent field variables for CaMCa
  CALL CMFE_Field_Initialise(CaMCaField,Err)
  CALL CMFE_EquationsSet_DependentCreateStart(CaMCaEquationsSet,CaMCaFieldUserNumber,CaMCaField,Err)
  CALL CMFE_Field_VariableLabelSet(CaMCaField,CMFE_FIELD_U_VARIABLE_TYPE,"CaMCa Field",Err)
  !Finish the equations set dependent field variables
  CALL CMFE_EquationsSet_DependentCreateFinish(CaMCaEquationsSet,Err)
  CALL CMFE_Field_ComponentValuesInitialise(CaMCaField,CMFE_FIELD_U_VARIABLE_TYPE, &
    & CMFE_FIELD_VALUES_SET_TYPE,1,init_CaMCa,Err)


  !Create the equations set material field variables - CaMCa
  !by default 2 comps for reac diff i.e. diff coeff in 1 direction set constant spatially = 1, and storage coeff set to 1
  CALL CMFE_Field_Initialise(CaMCaMaterialsField,Err)
  CALL CMFE_EquationsSet_MaterialsCreateStart(CaMCaEquationsSet,CaMCaMaterialsFieldUserNumber,CaMCaMaterialsField,Err)
  CALL CMFE_Field_VariableLabelSet(CaMCaMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,"CaMCa Materials Field",Err)
  !Finish the equations set materials field variables
  CALL CMFE_EquationsSet_MaterialsCreateFinish(CaMCaEquationsSet,Err)
  CALL CMFE_Field_ComponentValuesInitialise(CaMCaMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
   & 1,camcaDiffx,Err) !CaMCa diff coeff in x
  CALL CMFE_Field_ComponentValuesInitialise(CaMCaMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
   & 2,camcaDiffy,Err) !CaMCa diff coeff in y
  CALL CMFE_Field_ComponentValuesInitialise(CaMCaMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
   & 3,camcaDiffz,Err) !CaMCa diff coeff in z
  CALL CMFE_Field_ComponentValuesInitialise(CaMCaMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
   & 4,store_coeff,Err) ! storage coefficient

  CALL CMFE_Field_ParameterSetUpdateStart(CaMCaMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL CMFE_Field_ParameterSetUpdateFinish(CaMCaMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)

  !Set up source field for reaction diffusion equation set. Note that for the split problem subtype, the source field is not used at all.
  !iCaMCaField
  CALL CMFE_Field_Initialise(iCaMCaField,Err)
  CALL CMFE_EquationsSet_SourceCreateStart(CaMCaEquationsSet,iCaMCaFieldUserNumber,iCaMCaField,Err)
  CALL CMFE_Field_VariableLabelSet(iCaMCaField,CMFE_FIELD_U_VARIABLE_TYPE,"iCaMCa Field",Err)
  !Finish the equations set source field variables
  CALL CMFE_EquationsSet_SourceCreateFinish(CaMCaEquationsSet,Err)
  !Initialising the iCaField to zero everywhere. Modifying for RyRs in a later loop.
  CALL CMFE_Field_ComponentValuesInitialise(iCaMCaField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,0.0_CMISSDP,Err)

  !ATP equations
  CALL CMFE_EquationsSet_Initialise(ATPEquationsSet,Err)
  CALL CMFE_Field_Initialise(ATPEquationsSetField,Err)
  CALL CMFE_EquationsSet_CreateStart(ATPEquationsSetUserNumber,Region, & 
    & GeometricField,[CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMFE_EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE, &
    & CMFE_EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE], &
    & ATPEquationsSetFieldUserNumber,ATPEquationsSetField,ATPEquationsSet,Err)
  !Set the equations set to be a standard Diffusion no source problem
  !Finish creating the equations set
  CALL CMFE_EquationsSet_CreateFinish(ATPEquationsSet,Err)


  !Create the equations set dependent field variables for ATP
  CALL CMFE_Field_Initialise(ATPField,Err)
  CALL CMFE_EquationsSet_DependentCreateStart(ATPEquationsSet,ATPFieldUserNumber,ATPField,Err)
  CALL CMFE_Field_VariableLabelSet(ATPField,CMFE_FIELD_U_VARIABLE_TYPE,"ATP Field",Err)
  !Finish the equations set dependent field variables
  CALL CMFE_EquationsSet_DependentCreateFinish(ATPEquationsSet,Err)
  CALL CMFE_Field_ComponentValuesInitialise(ATPField,CMFE_FIELD_U_VARIABLE_TYPE, &
    & CMFE_FIELD_VALUES_SET_TYPE,1,init_ATP,Err)


  !Create the equations set material field variables - ATP
  !by default 2 comps for reac diff i.e. diff coeff in 1 direction set constant spatially = 1, and storage coeff set to 1
  CALL CMFE_Field_Initialise(ATPMaterialsField,Err)
  CALL CMFE_EquationsSet_MaterialsCreateStart(ATPEquationsSet,ATPMaterialsFieldUserNumber,ATPMaterialsField,Err)
  CALL CMFE_Field_VariableLabelSet(ATPMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,"ATP Materials Field",Err)
  !Finish the equations set materials field variables
  CALL CMFE_EquationsSet_MaterialsCreateFinish(ATPEquationsSet,Err)
  CALL CMFE_Field_ComponentValuesInitialise(ATPMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
   & 1,atpDiffx,Err) !ATP diff coeff in x
  CALL CMFE_Field_ComponentValuesInitialise(ATPMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
   & 2,atpDiffy,Err) !ATP diff coeff in y
  CALL CMFE_Field_ComponentValuesInitialise(ATPMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
   & 3,atpDiffz,Err) !ATP diff coeff in z
  CALL CMFE_Field_ComponentValuesInitialise(ATPMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
   & 4,store_coeff,Err) ! storage coefficient

  CALL CMFE_Field_ParameterSetUpdateStart(ATPMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL CMFE_Field_ParameterSetUpdateFinish(ATPMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)

  !Set up source field for reaction diffusion equation set. Note that for the split problem subtype, the source field is not used at all.
  !iATPField
  CALL CMFE_Field_Initialise(iATPField,Err)
  CALL CMFE_EquationsSet_SourceCreateStart(ATPEquationsSet,iATPFieldUserNumber,iATPField,Err)
  CALL CMFE_Field_VariableLabelSet(iATPField,CMFE_FIELD_U_VARIABLE_TYPE,"iATP Field",Err)
  !Finish the equations set source field variables
  CALL CMFE_EquationsSet_SourceCreateFinish(ATPEquationsSet,Err)
  !Initialising the iCaField to zero everywhere. Modifying for RyRs in a later loop.
  CALL CMFE_Field_ComponentValuesInitialise(iATPField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,0.0_CMISSDP,Err)

  !ATPCa equations
  CALL CMFE_EquationsSet_Initialise(ATPCaEquationsSet,Err)
  CALL CMFE_Field_Initialise(ATPCaEquationsSetField,Err)
  CALL CMFE_EquationsSet_CreateStart(ATPCaEquationsSetUserNumber,Region, & 
    & GeometricField,[CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMFE_EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE, &
    & CMFE_EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE], &
    & ATPCaEquationsSetFieldUserNumber,ATPCaEquationsSetField,ATPCaEquationsSet,Err)
  !Set the equations set to be a standard Diffusion no source problem
  !Finish creating the equations set
  CALL CMFE_EquationsSet_CreateFinish(ATPCaEquationsSet,Err)


  !Create the equations set dependent field variables for ATPCa
  CALL CMFE_Field_Initialise(ATPCaField,Err)
  CALL CMFE_EquationsSet_DependentCreateStart(ATPCaEquationsSet,ATPCaFieldUserNumber,ATPCaField,Err)
  CALL CMFE_Field_VariableLabelSet(ATPCaField,CMFE_FIELD_U_VARIABLE_TYPE,"ATPCa Field",Err)
  !Finish the equations set dependent field variables
  CALL CMFE_EquationsSet_DependentCreateFinish(ATPCaEquationsSet,Err)
  CALL CMFE_Field_ComponentValuesInitialise(ATPCaField,CMFE_FIELD_U_VARIABLE_TYPE, &
    & CMFE_FIELD_VALUES_SET_TYPE,1,init_ATPCa,Err)


  !Create the equations set material field variables - ATPCa
  !by default 2 comps for reac diff i.e. diff coeff in 1 direction set constant spatially = 1, and storage coeff set to 1
  CALL CMFE_Field_Initialise(ATPCaMaterialsField,Err)
  CALL CMFE_EquationsSet_MaterialsCreateStart(ATPCaEquationsSet,ATPCaMaterialsFieldUserNumber,ATPCaMaterialsField,Err)
  CALL CMFE_Field_VariableLabelSet(ATPCaMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,"ATPCa Materials Field",Err)
  !Finish the equations set materials field variables
  CALL CMFE_EquationsSet_MaterialsCreateFinish(ATPCaEquationsSet,Err)
  CALL CMFE_Field_ComponentValuesInitialise(ATPCaMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
   & 1,atpcaDiffx,Err) !ATPCa diff coeff in x
  CALL CMFE_Field_ComponentValuesInitialise(ATPCaMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
   & 2,atpcaDiffy,Err) !ATPCa diff coeff in y
  CALL CMFE_Field_ComponentValuesInitialise(ATPCaMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
   & 3,atpcaDiffz,Err) !ATPCa diff coeff in z
  CALL CMFE_Field_ComponentValuesInitialise(ATPCaMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
   & 4,store_coeff,Err) ! storage coefficient

  CALL CMFE_Field_ParameterSetUpdateStart(ATPCaMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL CMFE_Field_ParameterSetUpdateFinish(ATPCaMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)

  !Set up source field for reaction diffusion equation set. Note that for the split problem subtype, the source field is not used at all.
  !iATPCaField
  CALL CMFE_Field_Initialise(iATPCaField,Err)
  CALL CMFE_EquationsSet_SourceCreateStart(ATPCaEquationsSet,iATPCaFieldUserNumber,iATPCaField,Err)
  CALL CMFE_Field_VariableLabelSet(iATPCaField,CMFE_FIELD_U_VARIABLE_TYPE,"iATPCa Field",Err)
  !Finish the equations set source field variables
  CALL CMFE_EquationsSet_SourceCreateFinish(ATPCaEquationsSet,Err)
  !Initialising the iCaField to zero everywhere. Modifying for RyRs in a later loop.
  CALL CMFE_Field_ComponentValuesInitialise(iATPCaField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,0.0_CMISSDP,Err)
  CALL CMFE_Field_ParameterSetUpdateStart(iATPCaField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL CMFE_Field_ParameterSetUpdateFinish(iATPCaField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)

  !CaTnC
  CALL CMFE_Field_Initialise(CaTnCField,Err)
  CALL CMFE_Field_CreateStart(CaTnCFieldUserNumber,Region,CaTnCField,Err)
  CALL CMFE_Field_TypeSet(CaTnCField,CMFE_FIELD_GENERAL_TYPE,Err)
  CALL CMFE_Field_MeshDecompositionSet(CaTnCField,Decomposition,Err)
  CALL CMFE_Field_GeometricFieldSet(CaTnCField,GeometricField,Err)
  CALL CMFE_Field_NumberOfVariablesSet(CaTnCField,1,Err)
  CALL CMFE_Field_VariableTypesSet(CaTnCField,[CMFE_FIELD_U_VARIABLE_TYPE],Err)
  CALL CMFE_Field_DataTypeSet(CaTnCField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_DP_TYPE,Err)
  CALL CMFE_Field_DimensionSet(CaTnCField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_SCALAR_DIMENSION_TYPE,Err)
  CALL CMFE_Field_NumberOfComponentsSet(CaTnCField,CMFE_FIELD_U_VARIABLE_TYPE,1,Err)
  CALL CMFE_Field_VariableLabelSet(CaTnCField,CMFE_FIELD_U_VARIABLE_TYPE,"CaTnC Field",Err)
  CALL CMFE_Field_ComponentMeshComponentGet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE, & 
    & 1,GeometricMeshComponent,ERR)
  !Default to the geometric interpolation setup
  CALL CMFE_Field_ComponentMeshComponentSet(CaTnCField,CMFE_FIELD_U_VARIABLE_TYPE,1, &
    & GeometricMeshComponent,ERR)            
  !Specify the interpolation to be same as geometric interpolation
  CALL CMFE_Field_ComponentInterpolationSet(CaTnCField,CMFE_FIELD_U_VARIABLE_TYPE,1, &
    & CMFE_FIELD_NODE_BASED_INTERPOLATION,ERR)
  CALL CMFE_Field_CreateFinish(CaTnCField,Err)
  !Initialise CaTnC concentration to equilibrium value
  !Set the values to be nodally varying - mito nodes with different concentrations than myo regions
  CALL CMFE_Field_ComponentValuesInitialise(CaTnCField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,init_CaTnC,Err)
  CALL CMFE_Field_ParameterSetUpdateStart(CaTnCField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL CMFE_Field_ParameterSetUpdateFinish(CaTnCField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)


  !CaSRBuf
  CALL CMFE_Field_Initialise(CaSRBufField,Err)
  CALL CMFE_Field_CreateStart(CaSRBufFieldUserNumber,Region,CaSRBufField,Err)
  CALL CMFE_Field_TypeSet(CaSRBufField,CMFE_FIELD_GENERAL_TYPE,Err)
  CALL CMFE_Field_MeshDecompositionSet(CaSRBufField,Decomposition,Err)
  CALL CMFE_Field_GeometricFieldSet(CaSRBufField,GeometricField,Err)
  CALL CMFE_Field_NumberOfVariablesSet(CaSRBufField,1,Err)
  CALL CMFE_Field_VariableTypesSet(CaSRBufField,[CMFE_FIELD_U_VARIABLE_TYPE],Err)
  CALL CMFE_Field_DataTypeSet(CaSRBufField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_DP_TYPE,Err)
  CALL CMFE_Field_DimensionSet(CaSRBufField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_SCALAR_DIMENSION_TYPE,Err)
  CALL CMFE_Field_NumberOfComponentsSet(CaSRBufField,CMFE_FIELD_U_VARIABLE_TYPE,1,Err)
  CALL CMFE_Field_VariableLabelSet(CaSRBufField,CMFE_FIELD_U_VARIABLE_TYPE,"CaSRBuf Field",Err)
  CALL CMFE_Field_ComponentMeshComponentGet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE, & 
    & 1,GeometricMeshComponent,ERR)
  !Default to the geometric interpolation setup
  CALL CMFE_Field_ComponentMeshComponentSet(CaSRBufField,CMFE_FIELD_U_VARIABLE_TYPE,1, &
    & GeometricMeshComponent,ERR)            
  !Specify the interpolation to be same as geometric interpolation
  CALL CMFE_Field_ComponentInterpolationSet(CaSRBufField,CMFE_FIELD_U_VARIABLE_TYPE,1, &
    & CMFE_FIELD_NODE_BASED_INTERPOLATION,ERR)
  CALL CMFE_Field_CreateFinish(CaSRBufField,Err)
  !Initialise CaSRBuf concentration to equilibrium value
  !Set the values to be nodally varying - mito nodes with different concentrations than myo regions
  CALL CMFE_Field_ComponentValuesInitialise(CaSRBufField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,init_CaSRBuf,Err)
  CALL CMFE_Field_ParameterSetUpdateStart(CaSRBufField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL CMFE_Field_ParameterSetUpdateFinish(CaSRBufField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)

  !CaDyad
  CALL CMFE_Field_Initialise(CaDyadField,Err)
  CALL CMFE_Field_CreateStart(CaDyadFieldUserNumber,Region,CaDyadField,Err)
  CALL CMFE_Field_TypeSet(CaDyadField,CMFE_FIELD_GENERAL_TYPE,Err)
  CALL CMFE_Field_MeshDecompositionSet(CaDyadField,Decomposition,Err)
  CALL CMFE_Field_GeometricFieldSet(CaDyadField,GeometricField,Err)
  CALL CMFE_Field_NumberOfVariablesSet(CaDyadField,1,Err)
  CALL CMFE_Field_VariableTypesSet(CaDyadField,[CMFE_FIELD_U_VARIABLE_TYPE],Err)
  CALL CMFE_Field_DataTypeSet(CaDyadField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_DP_TYPE,Err)
  CALL CMFE_Field_DimensionSet(CaDyadField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_SCALAR_DIMENSION_TYPE,Err)
  CALL CMFE_Field_NumberOfComponentsSet(CaDyadField,CMFE_FIELD_U_VARIABLE_TYPE,1,Err)
  CALL CMFE_Field_VariableLabelSet(CaDyadField,CMFE_FIELD_U_VARIABLE_TYPE,"CaDyad Field",Err)
  CALL CMFE_Field_ComponentMeshComponentGet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE, & 
    & 1,GeometricMeshComponent,ERR)
  !Default to the geometric interpolation setup
  CALL CMFE_Field_ComponentMeshComponentSet(CaDyadField,CMFE_FIELD_U_VARIABLE_TYPE,1, &
    & GeometricMeshComponent,ERR)            
  !Specify the interpolation to be same as geometric interpolation
  CALL CMFE_Field_ComponentInterpolationSet(CaDyadField,CMFE_FIELD_U_VARIABLE_TYPE,1, &
    & CMFE_FIELD_NODE_BASED_INTERPOLATION,ERR)
  CALL CMFE_Field_CreateFinish(CaDyadField,Err)
  !Initialise CaDyad concentration to equilibrium value
  !Set the values to be nodally varying - mito nodes with different concentrations than myo regions
  CALL CMFE_Field_ComponentValuesInitialise(CaDyadField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,init_CaDyad,Err)
  CALL CMFE_Field_ParameterSetUpdateStart(CaDyadField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL CMFE_Field_ParameterSetUpdateFinish(CaDyadField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  
  !CaJSR
  CALL CMFE_Field_Initialise(CaJSRField,Err)
  CALL CMFE_Field_CreateStart(CaJSRFieldUserNumber,Region,CaJSRField,Err)
  CALL CMFE_Field_TypeSet(CaJSRField,CMFE_FIELD_GENERAL_TYPE,Err)
  CALL CMFE_Field_MeshDecompositionSet(CaJSRField,Decomposition,Err)
  CALL CMFE_Field_GeometricFieldSet(CaJSRField,GeometricField,Err)
  CALL CMFE_Field_NumberOfVariablesSet(CaJSRField,1,Err)
  CALL CMFE_Field_VariableTypesSet(CaJSRField,[CMFE_FIELD_U_VARIABLE_TYPE],Err)
  CALL CMFE_Field_DataTypeSet(CaJSRField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_DP_TYPE,Err)
  CALL CMFE_Field_DimensionSet(CaJSRField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_SCALAR_DIMENSION_TYPE,Err)
  CALL CMFE_Field_NumberOfComponentsSet(CaJSRField,CMFE_FIELD_U_VARIABLE_TYPE,1,Err)
  CALL CMFE_Field_VariableLabelSet(CaJSRField,CMFE_FIELD_U_VARIABLE_TYPE,"CaJSR Field",Err)
  CALL CMFE_Field_ComponentMeshComponentGet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE, & 
    & 1,GeometricMeshComponent,ERR)
  !Default to the geometric interpolation setup
  CALL CMFE_Field_ComponentMeshComponentSet(CaJSRField,CMFE_FIELD_U_VARIABLE_TYPE,1, &
    & GeometricMeshComponent,ERR)            
  !Specify the interpolation to be same as geometric interpolation
  CALL CMFE_Field_ComponentInterpolationSet(CaJSRField,CMFE_FIELD_U_VARIABLE_TYPE,1, &
    & CMFE_FIELD_NODE_BASED_INTERPOLATION,ERR)
  CALL CMFE_Field_CreateFinish(CaJSRField,Err)
  !Initialise CaJSR concentration to equilibrium value
  !Set the values to be nodally varying - mito nodes with different concentrations than myo regions
  CALL CMFE_Field_ComponentValuesInitialise(CaJSRField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,init_CaJSR,Err)
  CALL CMFE_Field_ParameterSetUpdateStart(CaJSRField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL CMFE_Field_ParameterSetUpdateFinish(CaJSRField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,Err)

  
!_____________________________________________________________________________________________________________________
  !Start to set up CellML Fields

  !Create the CellML environment
  CALL CMFE_CellML_Initialise(CellML,Err)
  CALL CMFE_CellML_CreateStart(CellMLUserNumber,Region,CellML,Err)
  !Import ryr release and buffer source model from a file
  CALL CMFE_CellML_ModelImport(CellML,RyRModel,ryrModelIndex,Err)
  ! set iCa as known so that it can be set as spatially varying in openCMFE_.
  !CALL CMFE_CellML_VariableSetAsKnown(CellML,ryrModelIndex,"CRU/iCa",Err)
  ! set RyRDensity as known so that it can be set as spatially varying in openCMFE_.
  CALL CMFE_CellML_VariableSetAsKnown(CellML,ryrModelIndex,"Dyad/NumRyR",Err)
  CALL CMFE_CellML_VariableSetAsKnown(CellML,ryrModelIndex,"Dyad/g_ilcc",Err)
  
  !to get from the CellML side. variables in cellml model that are not state variables, but are dependent on independent and state variables. 
  !- components of intermediate field
  !fluxes of the different buffers and CaRUs that I want to get out as intermediate variables
  CALL CMFE_CellML_VariableSetAsWanted(CellML,ryrModelIndex,"Dyad/J_ryr",Err)

  !Finish the CellML environment
  CALL CMFE_CellML_CreateFinish(CellML,Err)

    !Start the creation of CellML <--> OpenCMFE_ field maps
  !Mapping free calcium in openCMFE_ to that in cellml.

  CALL CMFE_CellML_FieldMapsCreateStart(CellML,Err)

  CALL CMFE_CellML_CreateFieldToCellMLMap(CellML,CaField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
    & ryrModelIndex,"Cytoplasm/Ca_cyto",CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL CMFE_CellML_CreateCellMLToFieldMap(CellML,ryrModelIndex,"Cytoplasm/Ca_cyto",CMFE_FIELD_VALUES_SET_TYPE, &
    & CaField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)

   !Mapping PopenField to P_open in the cellml model
  CALL CMFE_CellML_CreateFieldToCellMLMap(CellML,PopenField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
    & ryrModelIndex,"Dyad/P_open",CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL CMFE_CellML_CreateCellMLToFieldMap(CellML,ryrModelIndex,"Dyad/P_open",CMFE_FIELD_VALUES_SET_TYPE, &
    & PopenField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)

  CALL CMFE_CellML_CreateFieldToCellMLMap(CellML,iLCCField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
    & ryrModelIndex,"Dyad/g_ilcc",CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL CMFE_CellML_CreateCellMLToFieldMap(CellML,ryrModelIndex,"Dyad/g_ilcc",CMFE_FIELD_VALUES_SET_TYPE, &
    & iLCCField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)

   !Mapping NumRyRField to NumRyR in the cellml model
  CALL CMFE_CellML_CreateFieldToCellMLMap(CellML,NumRyRField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
    & ryrModelIndex,"Dyad/NumRyR",CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL CMFE_CellML_CreateCellMLToFieldMap(CellML,ryrModelIndex,"Dyad/NumRyR",CMFE_FIELD_VALUES_SET_TYPE, &
    & NumRyRField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)

  !Mapping Buffer-Complex resting values of cellml model to appropriate fields set up above

   !Mapping F
  CALL CMFE_CellML_CreateFieldToCellMLMap(CellML,FField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
    & ryrModelIndex,"Cytoplasm/Fluo_free",CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL CMFE_CellML_CreateCellMLToFieldMap(CellML,ryrModelIndex,"Cytoplasm/Fluo_free",CMFE_FIELD_VALUES_SET_TYPE, &
    & FField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)

   !Mapping FCa
  CALL CMFE_CellML_CreateFieldToCellMLMap(CellML,FCaField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
    & ryrModelIndex,"Cytoplasm/FluoCa",CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL CMFE_CellML_CreateCellMLToFieldMap(CellML,ryrModelIndex,"Cytoplasm/FluoCa",CMFE_FIELD_VALUES_SET_TYPE, &
    & FCaField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)

   !Mapping CaTnC
  CALL CMFE_CellML_CreateFieldToCellMLMap(CellML,CaTnCField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
    & ryrModelIndex,"Cytoplasm/CaTnC",CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL CMFE_CellML_CreateCellMLToFieldMap(CellML,ryrModelIndex,"Cytoplasm/CaTnC",CMFE_FIELD_VALUES_SET_TYPE, &
    & CaTnCField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)

  !Mapping CaM
  CALL CMFE_CellML_CreateFieldToCellMLMap(CellML,CaMField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
    & ryrModelIndex,"Cytoplasm/CaM_free",CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL CMFE_CellML_CreateCellMLToFieldMap(CellML,ryrModelIndex,"Cytoplasm/CaM_free",CMFE_FIELD_VALUES_SET_TYPE, &
    & CaMField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)

  !Mapping CaMCa
  CALL CMFE_CellML_CreateFieldToCellMLMap(CellML,CaMCaField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
    & ryrModelIndex,"Cytoplasm/CaMCa",CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL CMFE_CellML_CreateCellMLToFieldMap(CellML,ryrModelIndex,"Cytoplasm/CaMCa",CMFE_FIELD_VALUES_SET_TYPE, &
    & CaMCaField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)

  !Mapping ATP
  CALL CMFE_CellML_CreateFieldToCellMLMap(CellML,ATPField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
    & ryrModelIndex,"Cytoplasm/ATP_free",CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL CMFE_CellML_CreateCellMLToFieldMap(CellML,ryrModelIndex,"Cytoplasm/ATP_free",CMFE_FIELD_VALUES_SET_TYPE, &
    & ATPField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)

  !Mapping ATPCa
  CALL CMFE_CellML_CreateFieldToCellMLMap(CellML,ATPCaField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
    & ryrModelIndex,"Cytoplasm/ATPCa",CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL CMFE_CellML_CreateCellMLToFieldMap(CellML,ryrModelIndex,"Cytoplasm/ATPCa",CMFE_FIELD_VALUES_SET_TYPE, &
    & ATPCaField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)

  !Mapping Ca_dyad
  CALL CMFE_CellML_CreateFieldToCellMLMap(CellML,CaDyadField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
    & ryrModelIndex,"Dyad/Ca_dyad",CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL CMFE_CellML_CreateCellMLToFieldMap(CellML,ryrModelIndex,"Dyad/Ca_dyad",CMFE_FIELD_VALUES_SET_TYPE, &
    & CaDyadField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)

  !Mapping CaJSR
  CALL CMFE_CellML_CreateFieldToCellMLMap(CellML,CaJSRField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
    & ryrModelIndex,"SR/Ca_jsr",CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL CMFE_CellML_CreateCellMLToFieldMap(CellML,ryrModelIndex,"SR/Ca_jsr",CMFE_FIELD_VALUES_SET_TYPE, &
    & CaJSRField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)

  !Mapping CaSLBuf
  CALL CMFE_CellML_CreateFieldToCellMLMap(CellML,CaSLBufField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
    & ryrModelIndex,"Dyad/CaSLBuf",CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL CMFE_CellML_CreateCellMLToFieldMap(CellML,ryrModelIndex,"Dyad/CaSLBuf",CMFE_FIELD_VALUES_SET_TYPE, &
    & CaSLBufField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)

  !Mapping CaSRBuf
  CALL CMFE_CellML_CreateFieldToCellMLMap(CellML,CaSRBufField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
    & ryrModelIndex,"Cytoplasm/CaSRBuf",CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL CMFE_CellML_CreateCellMLToFieldMap(CellML,ryrModelIndex,"Cytoplasm/CaSRBuf",CMFE_FIELD_VALUES_SET_TYPE, &
    & CaSRBufField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)



  !Finish the creation of CellML <--> OpenCMFE_ field maps
  CALL CMFE_CellML_FieldMapsCreateFinish(CellML,Err)


  !Start the creation of the CellML models field. This field is an integer field that stores which nodes have which cellml model
  CALL CMFE_Field_Initialise(CellMLModelsField,Err)
  CALL CMFE_CellML_ModelsFieldCreateStart(CellML, &
    & CellMLModelsFieldUserNumber,CellMLModelsField,Err)
  !Finish the creation of the CellML models field
  CALL CMFE_CellML_ModelsFieldCreateFinish(CellML,Err)

  !By default all field parameters have default model value of 1, i.e. the first model. 
  ! assigning the bufferNryr cellml model (model 1) for all nodes.
  CALL CMFE_Field_ComponentValuesInitialise(CellMLModelsField, & 
    & CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,1_CMISSIntg,Err)
  !Start the creation of the CellML state field
  CALL CMFE_Field_Initialise(CellMLStateField,Err)
  CALL CMFE_CellML_StateFieldCreateStart(CellML, &
    & CellMLStateFieldUserNumber,CellMLStateField,Err)
  !Finish the creation of the CellML state field
  CALL CMFE_CellML_StateFieldCreateFinish(CellML,Err)

  !Start the creation of the CellML intermediate field
  CALL CMFE_Field_Initialise(CellMLIntermediateField,Err)
  CALL CMFE_CellML_IntermediateFieldCreateStart(CellML, &
    & CellMLIntermediateFieldUserNumber,CellMLIntermediateField,Err)
  !Finish the creation of the CellML intermediate field
  CALL CMFE_CellML_IntermediateFieldCreateFinish(CellML,Err)

  !Start the creation of CellML parameters field
  CALL CMFE_Field_Initialise(CellMLParametersField,Err)
  CALL CMFE_CellML_ParametersFieldCreateStart(CellML, &
    & CellMLParametersFieldUserNumber,CellMLParametersField,Err)
  !Finish the creation of CellML parameters
  CALL CMFE_CellML_ParametersFieldCreateFinish(CellML,Err)


  !Create the equations set equations for Ca
  WRITE(*,*) 'Creating the equations for the equation sets'
  CALL CMFE_Equations_Initialise(CaEquations,Err)
  CALL CMFE_EquationsSet_EquationsCreateStart(CaEquationsSet,CaEquations,Err)
  !Set the equations matrices sparsity type
  CALL CMFE_Equations_SparsityTypeSet(CaEquations,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  CALL CMFE_Equations_OutputTypeSet(CaEquations,CMFE_EQUATIONS_NO_OUTPUT,Err)
  !CALL CMFE_EquationsOutputTypeSet(Equations,CMFE_EquationsTimingOutput,Err)
  !CALL CMFE_EquationsOutputTypeSet(Equations,CMFE_EquationsMatrixOutput,Err)
  !CALL CMFE_EquationsOutputTypeSet(Equations,CMFE_EquationsElementMatrixOutput,Err)
  !Finish the equations set equations
  CALL CMFE_EquationsSet_EquationsCreateFinish(CaEquationsSet,Err)

  !Create the equations set equations for F
  CALL CMFE_Equations_Initialise(FEquations,Err)
  CALL CMFE_EquationsSet_EquationsCreateStart(FEquationsSet,FEquations,Err)
  !Set the equations matrices sparsity type
  CALL CMFE_Equations_SparsityTypeSet(FEquations,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  CALL CMFE_Equations_OutputTypeSet(FEquations,CMFE_EQUATIONS_NO_OUTPUT,Err)
  !CALL CMFE_EquationsOutputTypeSet(Equations,CMFE_EquationsTimingOutput,Err)
  !CALL CMFE_EquationsOutputTypeSet(Equations,CMFE_EquationsMatrixOutput,Err)
  !CALL CMFE_EquationsOutputTypeSet(Equations,CMFE_EquationsElementMatrixOutput,Err)
  !Finish the equations set equations
  CALL CMFE_EquationsSet_EquationsCreateFinish(FEquationsSet,Err)

  !Create the equations set equations for FCa
  CALL CMFE_Equations_Initialise(FCaEquations,Err)
  CALL CMFE_EquationsSet_EquationsCreateStart(FCaEquationsSet,FCaEquations,Err)
  !Set the equations matrices sparsity type
  CALL CMFE_Equations_SparsityTypeSet(FCaEquations,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  CALL CMFE_Equations_OutputTypeSet(FCaEquations,CMFE_EQUATIONS_NO_OUTPUT,Err)
  !CALL CMFE_EquationsOutputTypeSet(Equations,CMFE_EquationsTimingOutput,Err)
  !CALL CMFE_EquationsOutputTypeSet(Equations,CMFE_EquationsMatrixOutput,Err)
  !CALL CMFE_EquationsOutputTypeSet(Equations,CMFE_EquationsElementMatrixOutput,Err)
  !Finish the equations set equations
  CALL CMFE_EquationsSet_EquationsCreateFinish(FCaEquationsSet,Err)
  
  !Create the equations set equations for CaM
  CALL CMFE_Equations_Initialise(CaMEquations,Err)
  CALL CMFE_EquationsSet_EquationsCreateStart(CaMEquationsSet,CaMEquations,Err)
  !Set the equations matrices sparsity type
  CALL CMFE_Equations_SparsityTypeSet(CaMEquations,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  CALL CMFE_Equations_OutputTypeSet(CaMEquations,CMFE_EQUATIONS_NO_OUTPUT,Err)
  !CALL CMFE_EquationsOutputTypeSet(Equations,CMFE_EquationsTimingOutput,Err)
  !CALL CMFE_EquationsOutputTypeSet(Equations,CMFE_EquationsMatrixOutput,Err)
  !CALL CMFE_EquationsOutputTypeSet(Equations,CMFE_EquationsElementMatrixOutput,Err)
  !Finish the equations set equations
  CALL CMFE_EquationsSet_EquationsCreateFinish(CaMEquationsSet,Err)
  
  !Create the equations set equations for CaMCa
  CALL CMFE_Equations_Initialise(CaMCaEquations,Err)
  CALL CMFE_EquationsSet_EquationsCreateStart(CaMCaEquationsSet,CaMCaEquations,Err)
  !Set the equations matrices sparsity type
  CALL CMFE_Equations_SparsityTypeSet(CaMCaEquations,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  CALL CMFE_Equations_OutputTypeSet(CaMCaEquations,CMFE_EQUATIONS_NO_OUTPUT,Err)
  !CALL CMFE_EquationsOutputTypeSet(Equations,CMFE_EquationsTimingOutput,Err)
  !CALL CMFE_EquationsOutputTypeSet(Equations,CMFE_EquationsMatrixOutput,Err)
  !CALL CMFE_EquationsOutputTypeSet(Equations,CMFE_EquationsElementMatrixOutput,Err)
  !Finish the equations set equations
  CALL CMFE_EquationsSet_EquationsCreateFinish(CaMCaEquationsSet,Err)
  
  !Create the equations set equations for ATP
  CALL CMFE_Equations_Initialise(ATPEquations,Err)
  CALL CMFE_EquationsSet_EquationsCreateStart(ATPEquationsSet,ATPEquations,Err)
  !Set the equations matrices sparsity type
  CALL CMFE_Equations_SparsityTypeSet(ATPEquations,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  CALL CMFE_Equations_OutputTypeSet(ATPEquations,CMFE_EQUATIONS_NO_OUTPUT,Err)
  !CALL CMFE_EquationsOutputTypeSet(Equations,CMFE_EquationsTimingOutput,Err)
  !CALL CMFE_EquationsOutputTypeSet(Equations,CMFE_EquationsMatrixOutput,Err)
  !CALL CMFE_EquationsOutputTypeSet(Equations,CMFE_EquationsElementMatrixOutput,Err)
  !Finish the equations set equations
  CALL CMFE_EquationsSet_EquationsCreateFinish(ATPEquationsSet,Err)

  !Create the equations set equations for ATPCa
  CALL CMFE_Equations_Initialise(ATPCaEquations,Err)
  CALL CMFE_EquationsSet_EquationsCreateStart(ATPCaEquationsSet,ATPCaEquations,Err)
  !Set the equations matrices sparsity type
  CALL CMFE_Equations_SparsityTypeSet(ATPCaEquations,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  CALL CMFE_Equations_OutputTypeSet(ATPCaEquations,CMFE_EQUATIONS_NO_OUTPUT,Err)
  !CALL CMFE_EquationsOutputTypeSet(Equations,CMFE_EquationsTimingOutput,Err)
  !CALL CMFE_EquationsOutputTypeSet(Equations,CMFE_EquationsMatrixOutput,Err)
  !CALL CMFE_EquationsOutputTypeSet(Equations,CMFE_EquationsElementMatrixOutput,Err)
  !Finish the equations set equations
  CALL CMFE_EquationsSet_EquationsCreateFinish(ATPCaEquationsSet,Err)
  
  
!____________________________________________________________________________________________________________

  WRITE(*,*) 'Create the problem'
  CALL CMFE_Problem_Initialise(Problem,Err)
  CALL CMFE_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_CLASSICAL_FIELD_CLASS, &
    & CMFE_PROBLEM_REACTION_DIFFUSION_EQUATION_TYPE, &
    & CMFE_PROBLEM_CELLML_REAC_INTEG_REAC_DIFF_STRANG_SPLIT_SUBTYPE],Problem,Err)
  !Set the problem to be a strang split reaction diffusion problem
  !CALL CMFE_Problem_SpecificationSet(Problem,,Err)
  !Finish the creation of a problem.
  CALL CMFE_Problem_CreateFinish(Problem,Err)

  !Create the problem control
  CALL CMFE_Problem_ControlLoopCreateStart(Problem,Err)
  !Get the control loop
  CALL CMFE_ControlLoop_Initialise(ControlLoop,Err)
  CALL CMFE_Problem_ControlLoopGet(Problem,CMFE_CONTROL_LOOP_NODE,ControlLoop,Err)
  !Set the times
  CALL CMFE_ControlLoop_TimesSet(ControlLoop,startT,endT,Tstep,Err)
  CALL CMFE_ControlLoop_TimeOutputSet(ControlLoop,10,Err)
  CALL CMFE_ControlLoop_OutputTypeSet(ControlLoop,CMFE_CONTROL_LOOP_PROGRESS_OUTPUT,Err)
  !CALL CMFE_ControlLoopTimesSet(ControlLoop,0.0_CMISSDP,5.00_CMISSDP,0.01_CMISSDP,Err)
  !Finish creating the problem control loop
  CALL CMFE_Problem_ControlLoopCreateFinish(Problem,Err)

!______________________________________________________________________________________________________________
  WRITE(*,*) 'Set up the problem solvers for Strang splitting'
  !note, this example is contrived to have strang splitting, when it could be solved as a simple evaluation (as opposed to integration) of source and diffusion
  CALL CMFE_Problem_SolversCreateStart(Problem,Err)
  !First solver is a DAE solver
  CALL CMFE_Solver_Initialise(Solver,Err)
  CALL CMFE_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL CMFE_Solver_DAESolverTypeSet(Solver,CMFE_SOLVER_DAE_EULER,Err)
  CALL CMFE_Solver_LibraryTypeSet(Solver,CMFE_SOLVER_CMISS_LIBRARY,Err)
  CALL CMFE_Solver_DAETimeStepSet(Solver,ODE_TIME_STEP,Err)
  CALL CMFE_Solver_OutputTypeSet(Solver,CMFE_SOLVER_NO_OUTPUT,Err)

  !Second solver is the dynamic solver for solving the parabolic equation
  CALL CMFE_Solver_Initialise(Solver,Err)
  CALL CMFE_Solver_Initialise(LinearSolver,Err)
  CALL CMFE_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,2,Solver,Err)
  !set theta - backward vs forward time step parameter
  CALL CMFE_Solver_DynamicThetaSet(Solver,1.0_CMISSDP,Err)
  CALL CMFE_Solver_OutputTypeSet(Solver,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL CMFE_SolverOutputTypeSet(Solver,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
  !get the dynamic linear solver from the solver
  CALL CMFE_Solver_DynamicLinearSolverGet(Solver,LinearSolver,Err)
  !set linear solver to be direct solver. Note, I found this stuff in fluidmechanics/darcy/dynamic/src example
  !CALL CMFE_SolverLinearTypeSet(LinearSolver,CMFE_SolverLinearDirectSolveType,Err)
  !CALL CMFE_SolverLibraryTypeSet(LinearSolver,CMFE_SolverCMFE_Library,Err)
  !CALL CMFE_SolverLinearTypeSet(LinearSolver,CMFE_SolverLinearDirectSolveType,Err)
  !CALL CMFE_SolverLibraryTypeSet(LinearSolver,CMFE_SolverMUMPSLibrary,Err)
  CALL CMFE_Solver_LinearIterativeMaximumIterationsSet(LinearSolver,1000,Err)


  !Third solver is another DAE solver
  CALL CMFE_Solver_Initialise(Solver,Err)
  CALL CMFE_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,3,Solver,Err)
  CALL CMFE_Solver_DAESolverTypeSet(Solver,CMFE_SOLVER_DAE_EULER,Err)
  CALL CMFE_Solver_LibraryTypeSet(Solver,CMFE_SOLVER_CMISS_LIBRARY,Err)
  CALL CMFE_Solver_DAETimeStepSet(Solver,ODE_TIME_STEP,Err) !set the third solver's integration time step
  CALL CMFE_Solver_OutputTypeSet(Solver,CMFE_SOLVER_NO_OUTPUT,Err)
  !CALL CMFE_SolverOutputTypeSet(Solver,CMFE_SolverTimingOutput,Err)
  !CALL CMFE_SolverOutputTypeSet(Solver,CMFE_SolverSolverOutput,Err)
  !CALL CMFE_SolverOutputTypeSet(Solver,CMFE_SOLVER_PROGRESS_OUTPUT,Err)

  !Finish the creation of the problem solver
  CALL CMFE_Problem_SolversCreateFinish(Problem,Err)


  !Start the creation of the problem solver CellML equations
  CALL CMFE_Problem_CellMLEquationsCreateStart(Problem,Err)
  !Get the first solver  
  !Get the CellML equations
  CALL CMFE_Solver_Initialise(Solver,Err)
  CALL CMFE_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL CMFE_CellMLEquations_Initialise(CellMLEquations,Err)
  CALL CMFE_Solver_CellMLEquationsGet(Solver,CellMLEquations,Err)
  !Add in the CellML environement
  CALL CMFE_CellMLEquations_CellMLAdd(CellMLEquations,CellML,CellMLIndex,Err)

  !Get the third solver  
  !Get the CellML equations
  CALL CMFE_Solver_Initialise(Solver,Err)
  CALL CMFE_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,3,Solver,Err)
  CALL CMFE_CellMLEquations_Initialise(CellMLEquations,Err)
  CALL CMFE_Solver_CellMLEquationsGet(Solver,CellMLEquations,Err)
  !Add in the CellML environement
  CALL CMFE_CellMLEquations_CellMLAdd(CellMLEquations,CellML,CellMLIndex,Err)

  !Finish the creation of the problem solver CellML equations
  CALL CMFE_Problem_CellMLEquationsCreateFinish(Problem,Err)

  !Start the creation of the problem solver equations
  CALL CMFE_Problem_SolverEquationsCreateStart(Problem,Err)
  !Get the second solver  
  !Get the solver equations
  CALL CMFE_Solver_Initialise(Solver,Err)
  CALL CMFE_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,2,Solver,Err)
  CALL CMFE_SolverEquations_Initialise(SolverEquations,Err)
  CALL CMFE_Solver_SolverEquationsGet(Solver,SolverEquations,Err)
  !Set the solver equations sparsity
  !CALL CMFE_SolverEquationsSparsityTypeSet(SolverEquations,CMFE_SolverEquationsSparseMatrices,Err)
  CALL CMFE_SolverEquations_SparsityTypeSet(SolverEquations,CMFE_SOLVER_SPARSE_MATRICES,Err)  
  !Add in the equations set for Ca, F and FCa
  CALL CMFE_SolverEquations_EquationsSetAdd(SolverEquations,CaEquationsSet,EquationsSetIndex,Err)
  CALL CMFE_SolverEquations_EquationsSetAdd(SolverEquations,FEquationsSet,EquationsSetIndex,Err)
  CALL CMFE_SolverEquations_EquationsSetAdd(SolverEquations,FCaEquationsSet,EquationsSetIndex,Err)
  CALL CMFE_SolverEquations_EquationsSetAdd(SolverEquations,CaMEquationsSet,EquationsSetIndex,Err)
  CALL CMFE_SolverEquations_EquationsSetAdd(SolverEquations,CaMCaEquationsSet,EquationsSetIndex,Err)
  CALL CMFE_SolverEquations_EquationsSetAdd(SolverEquations,ATPEquationsSet,EquationsSetIndex,Err)
  CALL CMFE_SolverEquations_EquationsSetAdd(SolverEquations,ATPCaEquationsSet,EquationsSetIndex,Err)

  
  !Finish the creation of the problem solver equations
  CALL CMFE_Problem_SolverEquationsCreateFinish(Problem,Err)

!_________________________________________________________________________________________________________
  WRITE(*,*) 'Set up boundary conditions'
  CALL CMFE_BoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL CMFE_SolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)
  DO node=1,NUMBER_OF_NODES
    NODE_BD_LABEL = NodeNums(node,2)
    NODE_NUMBER = NodeNums(node,1)
    CALL CMFE_Decomposition_NodeDomainGet(Decomposition,NODE_NUMBER,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      IF(NODE_BD_LABEL.EQ.10) THEN !set no flux bc if node is a boundary node.
        CONDITION = CMFE_BOUNDARY_CONDITION_FIXED
        VALUE=0.0_CMISSDP
        CALL CMFE_BoundaryConditions_SetNode(BoundaryConditions,CaField, &
         & CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,CMFE_NO_GLOBAL_DERIV, &
         & NODE_NUMBER,1,CONDITION,VALUE,Err) !(neumann boundary condition - no flux)
        CALL CMFE_BoundaryConditions_SetNode(BoundaryConditions,FField, &
         & CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,CMFE_NO_GLOBAL_DERIV, &
         & NODE_NUMBER,1,CONDITION,VALUE,Err) !(neumann boundary condition - no flux)
        CALL CMFE_BoundaryConditions_SetNode(BoundaryConditions,FCaField, &
         & CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,CMFE_NO_GLOBAL_DERIV, &
         & NODE_NUMBER,1,CONDITION,VALUE,Err) !(neumann boundary condition - no flux)
        CALL CMFE_BoundaryConditions_SetNode(BoundaryConditions,CaMField, &
         & CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,CMFE_NO_GLOBAL_DERIV, &
         & NODE_NUMBER,1,CONDITION,VALUE,Err) !(neumann boundary condition - no flux)
        CALL CMFE_BoundaryConditions_SetNode(BoundaryConditions,CaMCaField, &
         & CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,CMFE_NO_GLOBAL_DERIV, &
         & NODE_NUMBER,1,CONDITION,VALUE,Err) !(neumann boundary condition - no flux)
        CALL CMFE_BoundaryConditions_SetNode(BoundaryConditions,ATPField, &
         & CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,CMFE_NO_GLOBAL_DERIV, &
         & NODE_NUMBER,1,CONDITION,VALUE,Err) !(neumann boundary condition - no flux)
        CALL CMFE_BoundaryConditions_SetNode(BoundaryConditions,ATPCaField, &
         & CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,CMFE_NO_GLOBAL_DERIV, &
         & NODE_NUMBER,1,CONDITION,VALUE,Err) !(neumann boundary condition - no flux)
         

      ENDIF
    ENDIF
  ENDDO
  CALL CMFE_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)
!__________________________________________________________________________________________________________
  !Solve the problem
  CALL CMFE_Problem_Solve(Problem,Err)

  IF(EXPORT_FIELD) THEN
    CALL CMFE_Fields_Initialise(Fields,Err)
    CALL CMFE_Fields_Create(Region,Fields,Err)
    CALL CMFE_Fields_NodesExport(Fields,"Ca_Cube_Solution","FORTRAN",Err)
    CALL CMFE_Fields_ElementsExport(Fields,"Ca_Cube_Solution","FORTRAN",Err)
    CALL CMFE_Fields_Finalise(Fields,Err)
  ENDIF 
  !Finialise CMFE_
  CALL CMFE_Finalise(Err)

  WRITE(*,'(A)') "Program successfully completed."
  
  STOP


END PROGRAM CUBE_SPARK

