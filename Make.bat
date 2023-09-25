:: !***************************************************************
:: ! This file is part of TREKIS-3
:: !***************************************************************

@echo off
setlocal EnableDelayedExpansion

:: read argument from user
   SET arg1=%1

::  in case of empty argument, assume no debug
   IF "%1"=="" (
      SET arg1=NODEBUG
   )
:: shorthand notations:
   IF "%1"=="db" (
      SET arg1=DEBUGOMP
   )
   IF "%1"=="DB" (
      SET arg1=DEBUGOMP
   )
   IF "%1"=="fast" (
      SET arg1=FAST
   )
   IF "%1"=="slow" (
      SET arg1=FAST
   )

   SET "Starline=************************************************************************************"
   echo %Starline%
   echo Started compilation: %date% %time%

:: List of all program files to be compiled
   SET "List_of_files=Universal_MC_for_SHI_MAIN.F90"

:: List compiler options and the name of the executable:
   IF /I %arg1%==DEBUGOMP (
      echo %Starline%
      echo Compiling with DEBUGOMP option, OpenMP but no optimizations are included
      echo Started at: %date% %time%
      echo %Starline%

      :: List compiler options
      SET "Compile_options=/F9999999999 /QxHost /QaxAVX  /fpp /Qopenmp /D OMP_inside /real-size:64 /debug:all /Od /check:all /traceback /gen-interfaces /warn:interfaces /check:bounds /fpe:0 /Qfp-stack-check /fp:precise /Qvec /standard-semantics"
      :: Set name of the executable:
      SET "Name_of_exe=TREKIS_DEBUG_OMP.exe"
    ) ELSE (
      IF /I %arg1%==FAST (
            echo %Starline%
            echo Compiling with FAST option, OpenMP, no optimizations, no debug
            echo Started at: %date% %time%
            echo %Starline%

            :: List compiler options
            SET "Compile_options=/F9999999999 /fpp /Qopenmp /D OMP_inside /Qmkl=parallel /real-size:64 /Od /fpe:0 /fp:precise /Qvec /standard-semantics"

            :: Set name of the executable:
            SET "Name_of_exe=TREKIS_OMP.exe"

            del *.pdb
      ) ELSE (
         echo %Starline%
         echo Compiling for release, OpenMP and optimizations are included
         echo Started at: %date% %time%
         echo %Starline%

         :: List compiler options
         SET "Compile_options=  /Qopenmp /D OMP_inside /O3 /fpp /Qvec /Qipo /real-size:64 /standard-semantics /F9999999999 "

         :: Set name of the executable:
         SET "Name_of_exe=TREKIS.exe"

         del *.pdb
      )
    )

:: Compile modules
   ifort.exe -c %Compile_options% %List_of_files%

   echo %Starline%
   echo Assembling the files into executable: %Name_of_exe%
   echo Started at: %date% %time%
   echo %Starline%

:: Assemble the code from all created obj-files
   ifort.exe %Compile_options% *.obj /exe:%Name_of_exe%

   echo %Starline%
::   echo Completed: %date% %time%
   echo The program %Name_of_exe% was created at %date% %time%
   echo %Starline%


:: Remove files that are no longer needed
del *.obj *.mod
