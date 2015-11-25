PRO PROGRESS, MSG, CUR, MAX 
  WRITEU, -1, STRING(FORMAT='(%"\R",A,": ",d5.1,"%")', $ 
              MSG, (FLOAT(CUR)/MAX*100.)) 
  IF CUR EQ MAX THEN PRINT, '' 
END 
