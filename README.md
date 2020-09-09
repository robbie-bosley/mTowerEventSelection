# mTowerEventSelection
Customisable software to perform event selection on mTower data from February 2020


# Quickly run an example of Event Selection
1. Make sure you have a working version of ROOT
2. Edit line 22 of 'Analyse_mTower_skeleton.sh' to your preferred Run Number. For this example, we will use Run 1336:
    Analyse_mTower(1336)
3. Simply run 'Analyse_mTower_skeleton.sh':

```bash    
./Analyse_mTower_skeleton.sh
```
You should see the various processors compile, including all the class descriptions in the classes/ folder, and the subprocessors in the subprocessors/ folder.
This compiling will likely throw some warnings about signed/unsigned int comparison, and some variables possibly being uninitialized - this is totally normal and nothing to worry about. 
After these processors have compiled, the code should then run. Expect it to take anywhere up to an hour, depending on how many events you're running over!

4. After the analysis code has finished running, you should find two files: results_Run_1336_C2_10p.root and Run_1336_eventsLeft.root.
Run_1336_eventsLeft.root contains data saved as a ROOT TTree.
results_Run_1336_C2_10p.root contains data saved as ROOT histograms.

5. You can check your output ROOT files using the ROOT TBrowser:

```bash
root -l results_Run_1336_C2_10p.root
new TBrowser
```

