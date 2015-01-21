
TST=dumb_write

### dumb write (test matrix)
$(TST)_new$(dbg): $(TST).c
	$(CC) $(DBG) $(PRF) $(NEW) $(WAR) $(STA) $(TST).c -o $(TST)_new$(dbg)

$(TST)_old$(dbg): $(TST).c
	$(CC) $(DBG) $(PRF) $(WAR) $(STA) $(TST).c -o $(TST)_old$(dbg)

$(TST)_new$(opt): $(TST).c
	$(CC) $(OPT) $(PRF) $(NEW) $(WAR) $(STA) $(TST).c -o $(TST)_new$(opt)

$(TST)_old$(opt): $(TST).c
	$(CC) $(OPT) $(PRF) $(WAR) $(STA) $(TST).c -o $(TST)_old$(opt)


