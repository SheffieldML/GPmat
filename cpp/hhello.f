	CALL S(12HHello, world, 12)
	END
	SUBROUTINE S(MSG,N)
	INTEGER K, N, M
	INTEGER MSG(1)
	M = (N + 3) / 4
	WRITE (6,'(20A4)') (MSG(K), K = 1,M)
	END
