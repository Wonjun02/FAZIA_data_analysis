0.find_PID -> Geant4 simulation에서 각 파티클의 dE-E 를 피팅, 저장
1.PID_dist -> find_PID로 계산한 피팅함수를 이용해 데이터의 PID index를 지정 + 가이드라인 시각화
2.fit_PID_dist -> Figure of Merit(FoM)을 계산하기 위해 pid_dist를 double gaussian + pol2로 피팅
