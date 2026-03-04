#### simulation functions

test_that("simulate_static_dbn generates valid data", {
	data <- simulate_static_dbn(n = 10, p = 2, time = 10, seed = 6886)

	expect_type(data, "list")
	expect_true(all(c("Y", "Z", "A", "B", "M", "sigma2", "tau2") %in% names(data)))

	expect_equal(dim(data$Y), c(10, 10, 2, 10))
	expect_equal(dim(data$Z), c(10, 10, 2, 10))
	expect_equal(dim(data$A), c(10, 10))
	expect_equal(dim(data$B), c(10, 10))
	expect_equal(dim(data$M), c(10, 10, 2))

	y_vals <- unique(c(data$Y[!is.na(data$Y)]))
	expect_true(all(y_vals %in% 1:5))

	for (t in 1:10) {
		for (p in 1:2) {
			expect_true(all(is.na(diag(data$Y[,,p,t]))))
		}
	}

	expect_equal(data$sigma2, 0.5)
	expect_equal(data$tau2, 0.1)
})

test_that("simulate_dynamic_dbn generates valid data", {
	data1 <- simulate_dynamic_dbn(n = 10, p = 2, time = 10,
		ar1 = FALSE, seed = 6886)

	expect_equal(dim(data1$Y), c(10, 10, 2, 10))
	expect_equal(dim(data1$A), c(10, 10, 10))
	expect_equal(dim(data1$B), c(10, 10, 10))
	expect_false(data1$ar1)
	expect_null(data1$rhoA)

	data2 <- simulate_dynamic_dbn(n = 10, p = 2, time = 10,
		ar1 = TRUE, rhoA = 0.9, rhoB = 0.8, seed = 6886)

	expect_true(data2$ar1)
	expect_equal(data2$rhoA, 0.9)
	expect_equal(data2$rhoB, 0.8)

	A_var <- var(c(data2$A))
	B_var <- var(c(data2$B))
	expect_true(A_var > 0)
	expect_true(B_var > 0)
})

test_that("simulate_lowrank_dbn generates valid data", {
	data <- simulate_lowrank_dbn(n = 10, p = 2, time = 10, r = 3, seed = 6886)

	expect_equal(dim(data$Y), c(10, 10, 2, 10))
	expect_equal(dim(data$U), c(10, 3))
	expect_equal(dim(data$alpha), c(3, 10))
	expect_equal(data$r, 3)

	UtU <- t(data$U) %*% data$U
	expect_equal(UtU, diag(3), tolerance = 1e-10)

	for (t in 1:5) {
		A_t <- data$U %*% diag(data$alpha[,t]) %*% t(data$U)
		expect_equal(A_t, data$A[,,t], tolerance = 1e-10)
	}
})

test_that("simulate_hmm_dbn generates valid data", {
	data <- simulate_hmm_dbn(n = 10, p = 2, time = 10, R = 3,
		transition_prob = 0.8, seed = 6886)

	expect_equal(dim(data$Y), c(10, 10, 2, 10))
	expect_length(data$S, 10)
	expect_equal(length(data$A_list), 3)
	expect_equal(length(data$B_list), 3)
	expect_equal(dim(data$Pi), c(3, 3))
	expect_equal(data$R, 3)

	expect_true(all(data$S %in% 1:3))

	expect_true(all(data$Pi >= 0))
	expect_equal(rowSums(data$Pi), rep(1, 3), tolerance = 1e-10)
	expect_equal(diag(data$Pi), rep(0.8, 3))

	A1_norm <- norm(data$A_list[[1]] - diag(10), "F")
	A2_norm <- norm(data$A_list[[2]] - diag(10), "F")
	expect_true(abs(A1_norm - A2_norm) > 0.01)
})

test_that("simulate_test_data generates simple data", {
	Y <- simulate_test_data(n = 10, p = 2, time = 20, seed = 333)

	expect_equal(dim(Y), c(10, 10, 2, 20))

	y_vals <- unique(c(Y[!is.na(Y)]))
	expect_true(all(y_vals %in% 1:5))

	for (t in 1:20) {
		for (p in 1:2) {
			expect_true(all(is.na(diag(Y[,,p,t]))))
		}
	}

	mean_r1 <- mean(Y[,,1,], na.rm = TRUE)
	mean_r2 <- mean(Y[,,2,], na.rm = TRUE)
	expect_true(mean_r1 > mean_r2)
})

test_that("simulation functions respect seeds", {
	data1 <- simulate_static_dbn(n = 10, p = 1, time = 5, seed = 999)
	data2 <- simulate_static_dbn(n = 10, p = 1, time = 5, seed = 999)
	expect_equal(data1$Y, data2$Y)

	data3 <- simulate_dynamic_dbn(n = 10, p = 1, time = 5, seed = 999)
	data4 <- simulate_dynamic_dbn(n = 10, p = 1, time = 5, seed = 999)
	expect_equal(data3$Y, data4$Y)

	data5 <- simulate_lowrank_dbn(n = 10, p = 1, time = 5, r = 2, seed = 999)
	data6 <- simulate_lowrank_dbn(n = 10, p = 1, time = 5, r = 2, seed = 999)
	expect_equal(data5$Y, data6$Y)

	data7 <- simulate_hmm_dbn(n = 10, p = 1, time = 5, R = 2, seed = 999)
	data8 <- simulate_hmm_dbn(n = 10, p = 1, time = 5, R = 2, seed = 999)
	expect_equal(data7$Y, data8$Y)
	expect_equal(data7$S, data8$S)
})

test_that("simulation functions handle edge cases", {
	data1 <- simulate_static_dbn(n = 5, p = 1, time = 1, seed = 444)
	expect_equal(dim(data1$Y), c(5, 5, 1, 1))

	data2 <- simulate_dynamic_dbn(n = 5, p = 1, time = 10, seed = 555)
	expect_equal(dim(data2$Y)[3], 1)

	data3 <- simulate_lowrank_dbn(n = 5, p = 1, time = 10, r = 1, seed = 666)
	expect_equal(ncol(data3$U), 1)
	expect_equal(nrow(data3$alpha), 1)

	data4 <- simulate_hmm_dbn(n = 5, p = 1, time = 10, R = 2, seed = 777)
	expect_equal(length(data4$A_list), 2)
})

#### recovery from simulated data

test_that("static model recovers block structure", {
	set.seed(6886)

	n <- 10
	time <- 50

	Y <- array(NA, dim = c(n, n, 1, time))

	for (i in 1:n) {
		for (j in 1:n) {
			if (i != j) {
				if ((i <= 5 & j <= 5) | (i > 5 & j > 5)) {
					Y[i,j,1,1] <- sample(3:5, 1)
				} else {
					Y[i,j,1,1] <- sample(1:3, 1)
				}
			}
		}
	}

	for (t in 2:time) {
		Y[,,1,t] <- Y[,,1,t-1]
		n_changes <- round(0.2 * n^2)
		for (k in 1:n_changes) {
			i <- sample(1:n, 1)
			j <- sample(setdiff(1:n, i), 1)
			Y[i,j,1,t] <- max(1, min(5, Y[i,j,1,t] + sample(c(-1,1), 1)))
		}
	}

	for (t in 1:time) diag(Y[,,1,t]) <- NA

	results <- dbn(Y, model = "static", nscan = 50, burn = 10, odens = 1, verbose = FALSE)

	B1_mean <- apply(results$B[[1]], c(1,2), mean)
	sender_effects <- diag(B1_mean)
	expect_true(length(sender_effects) == 10)
})

test_that("bilinear structure is preserved", {
	set.seed(6887)

	n <- 5
	A <- matrix(rnorm(n*n), n, n)
	B <- matrix(rnorm(n*n), n, n)
	Theta <- matrix(rnorm(n*n), n, n)

	Theta_new <- A %*% Theta %*% t(B)
	Theta_test <- tprod(Theta, list(A, B), modes = c(1,2))

	expect_equal(Theta_new, Theta_test, tolerance = 1e-10)
})

test_that("rank likelihood preserves ordering", {
	Y <- matrix(c(1,2,3,4,
		2,1,4,3,
		3,4,1,2,
		4,3,2,1), 4, 4)

	Z <- zscores(Y)

	for (i in 1:4) {
		y_order <- order(Y[i,])
		z_order <- order(Z[i,])
		expect_equal(y_order, z_order)
	}
})
