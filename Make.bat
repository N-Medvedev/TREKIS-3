:: !***************************************************************
:: ! This file is part of TREKIS-3
:: !***************************************************************

@echo off
setlocal EnableDelayedExpansion

:: Go into the directory with source files:
cd Source_files

:: read argument from user
   SET arg1=%1
   SET OMP=%2

::  in case of empty argument, assume no debug
   IF "%1"=="" (
      SET arg1=NODEBUG
   )
:: shorthand notations:
   IF "%1"=="debug" (
      SET arg1=DEBUGOMP
   )
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

   IF "%OMP%"=="" (
      SET OMP=YES
   )

   IF "%OMP%"=="NO" (
      SET OMP=no
   )

   SET "Starline=************************************************************************************"
   echo %Starline%
   echo Started compilation: %date% %time%

:: List of all program files to be compiled
   SET "List_of_files=Universal_MC_for_SHI_MAIN.f90"

:: List compiler options and the name of the executable:
   IF /I %arg1%==DEBUGOMP (
      echo %Starline%
      :: List compiler options
      :: Without OpenMP
      IF /I %OMP%==no (
         echo Compiling with DEBUG option, no OpenMP or optimizations are included
         echo Started at: %date% %time%
         SET "Compile_options=/F9999999999 /QxHost /QaxAVX /fpp /real-size:64 /debug:all /Od /Qipo /traceback /gen-interfaces /warn:interfaces /check:bounds /fpe:0 /fp:precise /Qvec /standard-semantics /assume:nofpe_summary"
         SET "Name_of_exe=TREKIS_DEBUG.exe"    :: Set name of the executable
      ) ELSE (    :: With OpenMP
         echo Compiling with DEBUG_OMP option, OpenMP but no optimizations are included
         echo Started at: %date% %time%
         SET "Compile_options=/F9999999999 /QxHost /QaxAVX /fpp /Qopenmp /real-size:64 /debug:all /O1 /Qipo /traceback /gen-interfaces /warn:interfaces /check:bounds /fpe:0 /fp:precise /Qvec /standard-semantics /assume:nofpe_summary"
         SET "Name_of_exe=TREKIS_DEBUG_OMP.exe"    :: Set name of the executable
      )
      echo %Starline%

    ) ELSE (
      IF /I %arg1%==FAST (
            echo %Starline%

            :: Without OpenMP
            IF /I %OMP%==no (
               echo Compiling with FAST option, but no OpenMP, optimizations, or debug
               echo Started at: %date% %time%
               echo %Starline%
               :: List compiler options
               SET "Compile_options=/F9999999999 /fpp /real-size:64 /O1 /fpe:0 /fp:fast /Qipo /Qopt-report /standard-semantics /assume:nofpe_summary"
               :: Set name of the executable:
               SET "Name_of_exe=TREKIS_no_OMP.exe"
            ) ELSE (    :: With OpenMP
               echo Compiling with FAST option, OpenMP, no optimizations, no debug
               echo Started at: %date% %time%
               echo %Starline%
               :: List compiler options
               SET "Compile_options=/F9999999999 /fpp /Qopenmp /real-size:64 /O1 /fpe:0 /fp:fast /Qipo /Qopt-report /standard-semantics /assume:nofpe_summary"
               :: Set name of the executable:
               SET "Name_of_exe=TREKIS_OMP.exe"
            )

            del *.pdb
      ) ELSE (
         echo %Starline%
         echo Compiling for release, OpenMP and optimizations are included
         echo Started at: %date% %time%
         echo %Starline%

         :: List compiler options
         SET "Compile_options=/F9999999999 /fpp /Qopenmp /real-size:64 /O3 /Qipo /standard-semantics /assume:nofpe_summary /static"

         :: Set name of the executable:
         SET "Name_of_exe=TREKIS.exe"

         del *.pdb
      )
    )

:: Compile modules
   ifx.exe -c %Compile_options% %List_of_files%

   echo %Starline%
   echo Assembling the files into executable: %Name_of_exe%
   echo Started at: %date% %time%
   echo %Starline%

:: Assemble the code from all created obj-files
   ifx.exe %Compile_options% *.obj /exe:%Name_of_exe%

   echo %Starline%
::   echo Completed: %date% %time%
   echo The program %Name_of_exe% was created at %date% %time%
   echo %Starline%

:: Remove files that are no longer needed
   del *.obj *.mod

:: Go back into the parent directory from the source files:
cd ..\

:: Copy exe file from the source directory into the parent directory:
xcopy source_files\%Name_of_exe% %Name_of_exe%* /Y /Q

:: Delete the exe file from the soure directory:
del source_files\%Name_of_exe%
