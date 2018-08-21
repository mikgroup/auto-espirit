rm log/nohup1.out log/nohup1.err
rm log/nohup2.out log/nohup2.err
rm log/nohup3.out log/nohup3.err

nohup nice matlab -nosplash -nodisplay -nojvm -r "experiment_1; exit;" > log/nohup1.out 2> log/nohup1.err < /dev/null &
nohup nice matlab -nosplash -nodisplay -nojvm -r "experiment_2; exit;" > log/nohup2.out 2> log/nohup2.err < /dev/null &
nohup nice matlab -nosplash -nodisplay -nojvm -r "experiment_3; exit;" > log/nohup3.out 2> log/nohup3.err < /dev/null &
