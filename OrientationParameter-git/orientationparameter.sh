#!/usr/bin/expect -f

set timeout -1
send -- "export GMX_MAXBACKUP=-1"

expect "$"
spawn mkdir GRO

set counter 0
set trajframes 10
set frame 25012

while { $counter < $trajframes } {
expect "$"
spawn gmx trjconv -f ../traj.xtc -s ../traj.tpr -o GRO/$counter.gro -b $frame -e $frame -pbc mol

expect "Select a group:"
send -- "0\r"

expect "$ "
spawn python3 orientationparameter.py
expect "Enter an index file name: "
send -- "pedot-heavy-all.ndx\r"
expect "Enter a gro file name: "
send -- "GRO/$counter.gro\r"


set counter [ expr $counter+1 ]
set frame [ expr $frame+104 ]
}

expect eof
