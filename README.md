Quickstart:

* setup root
```
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh

localSetupROOT 5.32.04-x86_64-slc5-gcc43-opt
```

* get this package
```
wget https://github.com/gerbaudo/WhSs2lTruth/archive/WhSs2lTruth-00-00-01.tar.gz
tar xzf WhSs2lTruth-00-00-01.tar.gz
mv WhSs2lTruth-WhSs2lTruth-00-00-01 WhSs2lTruth
```
* get SusyNtuple
```
wget https://github.com/gerbaudo/SusyNtuple/archive/SusyNtuple-00-01-16.tar.gz
tar xzf SusyNtuple-00-01-16.tar.gz
mv SusyNtuple-SusyNtuple-00-01-16 SusyNtuple
```
* get all the other dependencies
  Usually pulled in by `installMinimalSUSYTools.sh`, but here it's pre-pakaged.
  Note that this file is on gpatlas[1,2].
```
tar xzf /home/gerbaudo/tmp/WhSs2lTruth_dependencies_for_Ben.tgz
```

* Configure RootCore (needed only once)
```
cd RootCore
./configure
source scripts/setup.sh
cd ..
rc find_packages
```

* Compile everything (needed each time the code is modified)
```
rc build
```

* run
```
cd WhSs2lTruth/run/
run_truth_selector -i `cat filelist_notauhad_WH_2Lep_1.txt `
```

the output will be

```
...
Counts summary:
counter_input          : 99999
counter_event_cleaning : 56972
counter_EE             : 16904
counter_EM             : 25207
counter_MM             : 14897
counter_EE_SRSS1 : 429 counter_EE_SRSS2 : 207
counter_EM_SRSS1 : 608 counter_EM_SRSS2 : 283
counter_MM_SRSS1 : 433 counter_MM_SRSS2 : 160
```

davide.gerbaudo@gmail.com
Feb 2015