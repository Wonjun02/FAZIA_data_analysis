0.draw_PID를 이용해 proton line위에 피팅 지점을 표시하고 데이터 생성
1.cut_data로 실험데이터를 kaliveda의 beth_bloch formula로 피팅, 피팅결과와 histogram저장
2.fit_diff_y로 실험데이터를 fitting함수로 projection 후, gauissian fitting. 그 후 결과 저장
3.test_cust으로 1,2의 피팅 함수를 통해 raw data에서 proton 부분만 남도록 cut
4.draw_cut으로 결과 정리
