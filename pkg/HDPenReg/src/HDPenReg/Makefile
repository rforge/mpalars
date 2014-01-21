all:
	($(MAKE)  -f lars.mk)
	($(MAKE) -f lassoModels.mk)

clean:
	($(RM) lassoModels/*.o)
	($(RM) lars/*.o)
	(cd lib && $(RM) lib/*.a)
