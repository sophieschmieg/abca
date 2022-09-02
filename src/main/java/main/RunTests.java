package main;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import fields.finitefields.FiniteField;
import fields.finitefields.FiniteField.FFE;
import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.Ring.ExtendedEuclideanResult;
import fields.interfaces.Ring.FactorizationResult;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.local.CompleteDVRExtension;
import fields.local.FormalPowerSeries;
import fields.local.FormalPowerSeries.PowerSeries;
import fields.local.PAdicField;
import fields.local.PAdicField.PAdicNumber;
import fields.numberfields.NumberField;
import fields.numberfields.NumberField.NFE;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.CoordinateRing;
import fields.polynomials.CoordinateRing.CoordinateRingElement;
import fields.polynomials.Monomial;
import fields.polynomials.MultivariatePolynomialRing;
import fields.vectors.FreeModule;
import fields.vectors.Matrix;
import fields.vectors.MatrixAlgebra;
import fields.vectors.MatrixModule;
import fields.vectors.Vector;
import util.MiscAlgorithms;
import varieties.FunctionField;
import varieties.RationalFunction;
import varieties.curves.elliptic.EllipticCurve;
import varieties.projective.ProjectivePoint;

public class RunTests {

	public static void main(String[] args) {
		Rationals q = Rationals.q();
		Integers z = Integers.z();
		UnivariatePolynomialRing<IntE> integerRing = z.getUnivariatePolynomialRing();
		UnivariatePolynomial<IntE> t1 = integerRing.getPolynomial(z.getInteger(-5), z.getInteger(2), z.getInteger(8),
				z.getInteger(-3), z.getInteger(-3), z.getInteger(0), z.getInteger(1), z.getInteger(0), z.getInteger(1));
		UnivariatePolynomial<IntE> t2 = integerRing.getPolynomial(z.getInteger(21), z.getInteger(-9), z.getInteger(-4), z.getInteger(0),
				z.getInteger(5), z.getInteger(0), z.getInteger(3));
		UnivariatePolynomial<IntE> t3 = integerRing.getPolynomial(z.getInteger(0), z.getInteger(4), z.getInteger(4), z.getInteger(1));
		UnivariatePolynomial<IntE> t4 = integerRing.getPolynomial(z.getInteger(6), z.getInteger(5), z.getInteger(1));

		System.out.println(t1);
		System.out.println(integerRing.discriminant(t1));
		System.out.println(t2);
		System.out.println(integerRing.discriminant(t2));
		System.out.println(t3);
		System.out.println(integerRing.discriminant(t3));
		System.out.println(t4);
		System.out.println(integerRing.discriminant(t4));
//		UnivariatePolynomialRing.ExtendedResultantResult<IntE> result12 = integerRing.extendedResultant(t1, t2);
//		System.out.println("Resultant: " + result12.getResultant());
//		System.out.println("GCD:       " + result12.getGcd());
//		System.out.println("Coeff1:    " + result12.getCoeff1());
//		System.out.println("Coeff2:    " + result12.getCoeff2());
//		System.out.println(integerRing.add(integerRing.multiply(result12.getCoeff1(), t1), integerRing.multiply(result12.getCoeff2(), t2)));
//		UnivariatePolynomialRing.ExtendedResultantResult<IntE> result34 = integerRing.extendedResultant(t3, t4);
//		System.out.println("Resultant: " + result34.getResultant());
//		System.out.println("GCD:       " + result34.getGcd());
//		System.out.println("Coeff1:    " + result34.getCoeff1());
//		System.out.println("Coeff2:    " + result34.getCoeff2());
//		System.out.println(integerRing.add(integerRing.multiply(result34.getCoeff1(), t3), integerRing.multiply(result34.getCoeff2(), t4)));
		System.exit(0);
		MultivariatePolynomialRing<Fraction> qxy = (MultivariatePolynomialRing<Fraction>) AbstractPolynomialRing
				.getPolynomialRing(q, 2, Monomial.GREVLEX);
		Map<Monomial, Fraction> map1 = new TreeMap<>();
		map1.put(qxy.getMonomial(new int[] { 0, 0 }), q.one());
		map1.put(qxy.getMonomial(new int[] { 1, 0 }), q.one());
		map1.put(qxy.getMonomial(new int[] { 2, 0 }), q.one());
		Map<Monomial, Fraction> map2 = new TreeMap<>();
		map2.put(qxy.getMonomial(new int[] { 0, 0 }), q.one());
		map2.put(qxy.getMonomial(new int[] { 0, 1 }), q.one());
		map2.put(qxy.getMonomial(new int[] { 0, 3 }), q.one());
		List<Polynomial<Fraction>> plist = new ArrayList<>();
		plist.add(qxy.getPolynomial(map1));
		plist.add(qxy.getPolynomial(map2));
		CoordinateRing<Fraction> cr1 =  qxy.getIdeal(plist).divideOut();
		CoordinateRingElement<Fraction> x2y = cr1.add(cr1.getEmbedding(qxy.getVar(1)),
				cr1.multiply(2, cr1.getEmbedding(qxy.getVar(2))));
		CoordinateRingElement<Fraction> power = cr1.one();
		List<Vector<Fraction>> asVectors = new ArrayList<>();
		for (int i = 0; i <= 6; i++) {
			System.out.println(power);
			asVectors.add(qxy.asVector(power.getElement(), new int[] { 2, 3 }));
			power = cr1.multiply(power, x2y);
		}
		Matrix<Fraction> asMatrix = Matrix.fromColumns(asVectors);
		System.out.println(asMatrix.toString());
		MatrixModule<Fraction> mm = new MatrixModule<>(q, 6, 7);
		System.out.println(mm.kernelBasis(asMatrix));
		UnivariatePolynomial<Fraction> mipo = q.getUnivariatePolynomialRing()
				.getPolynomial(mm.kernelBasis(asMatrix).get(0).asList());
		System.out.println(mipo);
		asVectors.remove(asVectors.size() - 1);
		Vector<Fraction> vectorX = qxy.asVector(qxy.getVar(1), new int[] { 2, 3 });
		Matrix<Fraction> solveXMatrix = Matrix.fromColumns(asVectors);
		System.out.println(solveXMatrix.toString());
		MatrixModule<Fraction> mm6 = new MatrixModule<>(q, 6, 6);
		NumberField numberfield = NumberField.getNumberField(mipo);
		System.out.println(numberfield);
		UnivariatePolynomialRing<NFE> nfring = numberfield.getUnivariatePolynomialRing();
		UnivariatePolynomial<NFE> xAsPoly = nfring.getEmbedding(
				q.getUnivariatePolynomialRing().getPolynomial(mm6.solve(solveXMatrix, vectorX).asList()),
				numberfield.getEmbeddingMap());
		Polynomial<NFE> miponfe = nfring.getEmbedding(mipo, numberfield.getEmbeddingMap());
		NFE alpha = nfring.evaluate(xAsPoly, numberfield.alpha());
		System.out.println(alpha);
		System.out.println(numberfield.add(numberfield.multiply(alpha, alpha), alpha, numberfield.one()));
		System.out.println(nfring.normalize(nfring.substitute(miponfe,
				Collections.singletonList(nfring.getPolynomial(alpha, numberfield.getInteger(2))))));
		System.out.println(numberfield
				.factorization(nfring.getPolynomial(numberfield.one(), numberfield.one(), numberfield.one())));
		PAdicField z3 = new PAdicField(BigInteger.valueOf(3), 4);
		UnivariatePolynomialRing<PAdicNumber> z3ring = z3.getUnivariatePolynomialRing();
//		 UnivariatePolynomial<PAdicBaseNumber> t =
//		 z3ring.getPolynomial(z3.getInteger(9), z3.getInteger(3),
//		 z3.getInteger(3), z3.getInteger(1));
		UnivariatePolynomial<PAdicNumber> t = z3ring.getPolynomial(z3.getInteger(9), z3.getInteger(9), z3.getInteger(6),
				z3.getInteger(3), z3.getInteger(1));
//		UnivariatePolynomial<PAdicBaseNumber> t = z3ring.getPolynomial(z3.getInteger(-9), z3.getInteger(0),
//				z3.getInteger(1));
//		UnivariatePolynomial<PAdicBaseNumber> t = z3ring.getPolynomial(z3.getInteger(3), z3.getInteger(3),
//				z3.getInteger(3), z3.getInteger(1));
	FactorizationResult<Polynomial<PAdicNumber>, PAdicNumber> factorss = z3.factorization(t);
		System.out.println(t);
		for (Polynomial<PAdicNumber> f : factorss.primeFactors()) {
			System.out.println("Factor: " + f);
		}
		UnivariatePolynomial<PAdicNumber> lifted = z3ring.getVar();
		System.out.println(t);
		System.out.println(lifted);
		CoordinateRing<PAdicNumber> cr = z3ring.getIdeal(Collections.singletonList(t)).divideOut();
		ExtendedEuclideanResult<Polynomial<PAdicNumber>> ee = z3ring.extendedEuclidean(lifted, t);
		if (!ee.getGcd().equals(z3ring.one())) {
			throw new ArithmeticException("Trivial factors not removed or extended Euclidean not normalized!");
		}
		CoordinateRingElement<PAdicNumber> invertX = cr.getEmbedding(ee.getCoeff1());
		invertX = cr.multiply(invertX, invertX);
		System.out.println(invertX);
		invertX = cr.multiply(3, invertX);
		// invertX = cr.add(invertX, cr.one());
		for (int i = 0; i < 20; i++) {
			System.out.println("Power: " + cr.power(invertX, i));
			System.out.println("gcd:   " + z3ring.gcd(t, cr.power(invertX, i).getElement()));
		}
		System.exit(0);

		PrimeField f5 = PrimeField.getPrimeField(5);
		FormalPowerSeries<PFE> powerSeries = new FormalPowerSeries<>(f5, 10);
		System.out.println(powerSeries.zero());
		System.out.println(powerSeries.one());
		System.out.println(powerSeries.uniformizer());
		System.out.println(powerSeries.inverse(powerSeries.uniformizer()));
		PowerSeries<PFE> xmone = powerSeries
				.getEmbedding(f5.getUnivariatePolynomialRing().getPolynomial(f5.negative(f5.one()), f5.one()));
		System.out.println(xmone);
		System.out.println(powerSeries.inverse(xmone));
		System.out.println(powerSeries.multiply(xmone, powerSeries.inverse(xmone)));
		System.out.println(powerSeries.add(xmone, powerSeries.one()));
		System.out.println(powerSeries.add(xmone, powerSeries.uniformizer()));
		System.out.println(powerSeries.add(xmone, powerSeries.inverse(powerSeries.uniformizer())));
		System.out.println();

		System.out.println(powerSeries.roundToRationalFunction(powerSeries.zero()));
		System.out.println(powerSeries.roundToRationalFunction(powerSeries.one()));
		System.out.println(powerSeries.roundToRationalFunction(powerSeries.uniformizer()));
		System.out.println(powerSeries.roundToRationalFunction(xmone));
		System.out.println(powerSeries.roundToRationalFunction(powerSeries.inverse(powerSeries.uniformizer())));
		System.out.println(powerSeries.roundToRationalFunction(powerSeries.inverse(xmone)));
		System.out.println(powerSeries.roundToRationalFunction(powerSeries.divide(powerSeries.uniformizer(), xmone)));
		PolynomialRing<PFE> p = AbstractPolynomialRing.getPolynomialRing(f5, 2, Monomial.GREVLEX);
		// Polynomial<PFE> t = p.subtract(p.getVarPower(1, 2), p.getVarPower(2, 2));
		// t = p.add(t, p.multiply(2, p.getVar(2)), p.getEmbedding(f5.getElement(-1)));
		// List<Polynomial<PFE>> fs = p.factorization(t);
		// for (Polynomial<PFE> f : fs) {
		// System.out.println(f);
		// }
//		Polynomial<PFE> t2 = p.subtract(
//				p.multiply(p.getVar(1), p.add(p.getVar(1), p.one()), p.subtract(p.getVar(1), p.one())),
//				p.getVarPower(2, 2));
//		List<Polynomial<PFE>> fs2 = p.factorization(t2);
//		for (Polynomial<PFE> f : fs2) {
//			System.out.println(f);
//		}
//		System.exit(0);
		// new ProduceSivOutput();
		PAdicField z5 = new PAdicField(BigInteger.valueOf(5), 12);
		// Rationals q = Rationals.q();
//		for (int i = 0; i < 100; i++) {
//			for (int j = 0; j < 100; j++) {
//				PAdicNumber n = z5.getInteger(i);
//				PAdicNumber m = z5.getInteger(j);
//				if (z5.toRational(n).getNumerator().getValue().intValueExact() != i) {
//					System.out.println(n + " != " + i);
//				}
//				if (z5.toRational(m).getNumerator().getValue().intValueExact() != j) {
//					System.out.println(n + " != " + i);
//				}
//				if (z5.toRational(z5.add(n, m)).getNumerator().getValue().intValueExact() != i + j) {
//					System.out.println(n + " + " + m + " != " + i + " + " + j);
//				}
//				if (z5.toRational(z5.multiply(n, m)).getNumerator().getValue().intValueExact() != i * j) {
//					System.out.println(n + " * " + m + " = " + z5.multiply(n, m) + " != " + i + " * " + j);
//				}
//				if (j != 0 && !z5.multiply(m, z5.divide(n, m)).equals(n)) {
//					System.out.println(n + " * (" + m + "/" + n + ") = " + n + " * " + z5.divide(m, n) + " = " + z5.multiply(m, z5.divide(n, m)) + " != " + n);
//				}
//			}
//		}
//		System.out.println();
//		PolynomialRing<Fraction> q4 = AbstractPolynomialRing.getPolynomialRing(q, 4, Monomial.LEX);
//		Polynomial<Fraction> p1 = q4.getVarPower(1, 2);
//		p1 = q4.subtract(p1, q4.getVarPower(2, 2));
//		p1 = q4.subtract(p1, q4.getVarPower(3, 2));
//		p1 = q4.add(p1, q4.getVarPower(4, 2));
//		p1 = q4.subtract(p1, q4.getInteger(3));
//		Polynomial<Fraction> p2 = q4.negative(q4.getVarPower(2, 2));
//		p2 = q4.add(p2, q4.getVarPower(4, 2));
//		p2 = q4.add(p2, q4.multiply(2, q4.getVar(1), q4.getVar(2)));
//		p2 = q4.subtract(p2, q4.multiply(2, q4.getVar(3), q4.getVar(4)));
//		Polynomial<Fraction> p3 = q4.multiply(2, q4.getVar(1), q4.getVar(3));
//		p3 = q4.subtract(p3, q4.multiply(2, q4.getVar(2), q4.getVar(4)));
//		Polynomial<Fraction> p4 = q4.multiply(2, q4.getVar(1), q4.getVar(4));
//		p4 = q4.add(p4, q4.multiply(2, q4.getVar(2), q4.getVar(3)));
//		p4 = q4.subtract(p4, q4.multiply(2, q4.getVar(2), q4.getVar(4)));
//		System.out.println(p1);
//		System.out.println(p2);
//		System.out.println(p3);
//		System.out.println(p4);
//		List<Polynomial<Fraction>> list = new ArrayList<>();
//		list.add(p1);
//		list.add(p2);
//		list.add(p3);
//		list.add(p4);
//		List<AffinePoint<Fraction>> solution = q4.solve(list);
//		System.out.println("Found " + solution.size() + " solutions");
//		for (AffinePoint<Fraction> point : solution) {
//			System.out.println(point);
//		}
//		Polynomial<Fraction> test = q4.subtract(q4.getVarPower(1, 2), q4.getVarPower(2, 2));
//		System.out.println(test);
//		List<Polynomial<Fraction>> testFactors = q4.factorization(test);
//		for (Polynomial<Fraction> factor : testFactors) {
//			System.out.println(factor);
//		}
//		
//		Ideal<Polynomial<Fraction>> idealOfInterest = q4.getIdeal(list);
//		System.out.println();
//		for (Polynomial<Fraction> polynomial : idealOfInterest.generators()) {
//			System.out.println(polynomial);
//		}
//		System.exit(0);
//		BigInteger order = BigInteger.ONE.shiftLeft(128).subtract(BigInteger.ONE).divide(BigInteger.ONE.shiftLeft(8).subtract(BigInteger.ONE));
//		System.out.println(order);
//		BigInteger coorder = BigInteger.ONE.shiftLeft(8).subtract(BigInteger.ONE);
//		System.out.println(coorder.toString(16));
//		System.out.println(f2to128.power(f2to128.alpha(), BigInteger.ONE.shiftLeft(128)));
//		System.out.println(f2to128.power(f2to128.alpha(), order));
//		FFE possible = f2to128.power(f2to128.alpha(), coorder);
//		System.out.println(possible);
//		System.out.println(f2to128.minimalPolynomial(possible));
//		List<Vector<PrimeFieldElement>> vectors = new ArrayList<>();
//		for (int i = 0; i <= 16; i++) {
//		vectors.add(f2to128.asVector(f2to128.power(possible, i)));
//		}
//		MatrixModule<PrimeFieldElement> module = new MatrixModule<>(f2, 128, 17);
//		List<Vector<PrimeFieldElement>> kernelBasis = module.kernelBasis(Matrix.fromColumns(vectors));
//		for (Vector<PrimeFieldElement> b : kernelBasis) {
//			System.out.println(b);
//		}
		// System.out.println(f2to128.asMatrix(possible));
		PAdicField fiveadic = new PAdicField(BigInteger.valueOf(5), 5);
		PAdicNumber five = fiveadic.getInteger(5);
		PAdicNumber four = fiveadic.getInteger(4);
		PAdicNumber minusone = fiveadic.getInteger(-1);
		PAdicNumber onefifth = fiveadic.inverse(fiveadic.getInteger(5));
		PAdicNumber onequarter = fiveadic.inverse(fiveadic.getInteger(4));
		PAdicNumber threequarter = fiveadic.divide(fiveadic.getInteger(3), fiveadic.getInteger(4));
		System.out.println("5 = " + five + " = " + fiveadic.toRational(five));
		System.out.println("4 = " + four + " = " + fiveadic.toRational(four));
		System.out.println("-1 = " + minusone + " = " + fiveadic.toRational(minusone));
		System.out.println("1/5 = " + onefifth + " = " + fiveadic.toRational(onefifth));
		System.out.println("1/4 = " + onequarter + " = " + fiveadic.toRational(onequarter));
		System.out.println("3/4 = " + threequarter + " = " + fiveadic.toRational(threequarter));

//		CompleteDVRExtension<PAdicNumber, PFE, PFE, FFE, FiniteField> z5roots = fiveadic
//				.getExtension(fiveadic.getUnivariatePolynomialRing().getPolynomial(fiveadic.one(), fiveadic.one(),
//						fiveadic.one(), fiveadic.one(), fiveadic.one()))
//				.extension();
//		System.out.println(z5roots.ramificationIndex());
//		System.out.println(z5roots.residueDegree());
//		System.out.println(z5roots.minimalPolynomial());
		/*
		 * System.out.println(z5roots.ramificationPolynomial());
		 * System.out.println(z5roots.inertiaPolynomial());
		 * System.out.println(z5roots.ramificationAlpha());
		 * System.out.println(z5roots.inertiaAlpha());
		 */
		UnivariatePolynomialRing<Fraction> qr = q.getUnivariatePolynomialRing();
		Polynomial<Fraction> poly = qr.getPolynomial(q.one(), q.zero(), q.zero(), q.zero(), q.zero(), q.zero(),
				q.zero(), q.one(), q.one());
		FactorizationResult<Polynomial<Fraction>, Polynomial<Fraction>> factors = qr.uniqueFactorization(poly);
		System.out.println("Factorization over the fractions:\n" + factors);
	
		Reals r = Reals.r(1024);
		Polynomial<Real> polyReal = MiscAlgorithms.mapPolynomial(poly, new MathMap<>() {

			@Override
			public Real evaluate(Fraction t) {
				return r.getEmbedding(t);
			}
		}, r.getUnivariatePolynomialRing());
		FactorizationResult<Polynomial<Real>, Polynomial<Real>> factorsReal = r.getUnivariatePolynomialRing().uniqueFactorization(polyReal);
		System.out.println("Factorization over the reals:\n" + factorsReal);
	
		UnivariatePolynomial<Fraction> minpoly = qr.getPolynomial(q.one(), q.one(), q.one());
		NumberField nf = NumberField.getNumberField(minpoly);
		Polynomial<NFE> polyNF = MiscAlgorithms.mapPolynomial(poly, new MathMap<>() {

			@Override
			public NFE evaluate(Fraction t) {
				return nf.getEmbedding(t);
			}
		}, nf.getUnivariatePolynomialRing());
		FactorizationResult<Polynomial<NFE>, Polynomial<NFE>> factorsNF = nf.getUnivariatePolynomialRing().uniqueFactorization(polyNF);
		System.out.println("Factorization over NumberField:\n" + factorsNF);
	
//		Map<Monomial, Fraction> intmap = new TreeMap<>(qr.getComparator());
//		intmap.put(qr.getMonomial(new int[] { 128 }), q.one());
//		intmap.put(qr.getMonomial(new int[] { 112 }), q.negative(q.one()));
//		intmap.put(qr.getMonomial(new int[] { 80 }), q.one());
//		intmap.put(qr.getMonomial(new int[] { 64 }), q.negative(q.one()));
//		intmap.put(qr.getMonomial(new int[] { 48 }), q.one());
//		intmap.put(qr.getMonomial(new int[] { 16 }), q.negative(q.one()));
//		intmap.put(qr.getMonomial(new int[] { 0 }), q.one());
//		Polynomial<Fraction> poly2 = qr.getPolynomial(intmap);
//		System.out.println("Factorization over Rational:");
//		List<Polynomial<Fraction>> factors2 = q.factorization(poly2);
//		for (Polynomial<Fraction> factor : factors2) {
//			System.out.println("factor: " + factor);
//		}

		NumberField extended = NumberField.getNumberField(minpoly);
		MatrixAlgebra<Fraction> algebra = new FreeModule<>(q, 4).matrixAlgebra();
		List<List<Fraction>> m = new ArrayList<>();
		List<Fraction> row1 = new ArrayList<>();
		row1.add(q.getInteger(1));
		row1.add(q.getInteger(-1));
		row1.add(q.getInteger(11));
		row1.add(q.getInteger(-71));
		m.add(row1);
		List<Fraction> row2 = new ArrayList<>();
		row2.add(q.getInteger(1));
		row2.add(q.getInteger(-3));
		row2.add(q.getInteger(7));
		row2.add(q.getInteger(-49));
		m.add(row2);
		List<Fraction> row3 = new ArrayList<>();
		row3.add(q.getInteger(0));
		row3.add(q.getInteger(3));
		row3.add(q.getInteger(-10));
		row3.add(q.getInteger(35));
		m.add(row3);
		List<Fraction> row4 = new ArrayList<>();
		row4.add(q.getInteger(0));
		row4.add(q.getInteger(1));
		row4.add(q.getInteger(-10));
		row4.add(q.getInteger(37));
		m.add(row4);
		Matrix<Fraction> matrix = new Matrix<>(m);
		System.out.println(matrix);
		List<Fraction> rhs = new ArrayList<>();
		rhs.add(q.getInteger(1));
		rhs.add(q.getInteger(0));
		rhs.add(q.getInteger(0));
		rhs.add(q.getInteger(0));
		Vector<Fraction> vector = new Vector<>(rhs);
		System.out.println(vector);
//		Vector<Fraction> result = algebra.solve(matrix, vector);
//		System.out.println(result);
		alpha = extended.alpha();
		NFE iE = extended.multiply(extended.getEmbedding(q.getFraction(z.getInteger(1), z.getInteger(24))),
				extended.add(extended.multiply(-127, alpha),
						extended.add(extended.multiply(-5, extended.power(alpha, 3)),
								extended.multiply(-19, extended.power(alpha, 5)),
								extended.multiply(-5, extended.power(alpha, 7)))));
		NFE fourthRootOfTwo = extended.subtract(alpha, iE);
		NFE sqrtOfTwo = extended.multiply(fourthRootOfTwo, fourthRootOfTwo);
		System.out.println(iE);
		System.out.println(fourthRootOfTwo);
		System.out.println(extended.multiply(iE, iE));
		System.out.println(extended.power(fourthRootOfTwo, 4));
		System.out.println("N((i-1)/sqrt(2)) = "
				+ extended.norm(extended.divide(extended.subtract(iE, extended.one()), sqrtOfTwo)));
		System.out.println("N(i) = " + extended.norm(iE));
		System.out.println("N(i-1) = " + extended.norm(extended.subtract(iE, extended.one())));
		System.out.println("N(4sqrt(2)) = " + extended.norm(fourthRootOfTwo));
		System.out.println("N(sqrt(2)) = " + extended.norm(sqrtOfTwo));
		System.out.println("N(2) = " + extended.norm(extended.getInteger(2)));
		System.out
				.println(
						"N((i-1-sqrt(2)-2) = "
								+ extended
										.norm(extended.divide(
												extended.subtract(iE,
														extended.add(extended.one(), sqrtOfTwo,
																extended.power(fourthRootOfTwo, 3))),
												extended.getInteger(2))));

		// System.exit(0);

		PrimeField base = PrimeField.getPrimeField(5);
		UnivariatePolynomialRing<PFE> ring = base.getUnivariatePolynomialRing();
		UnivariatePolynomial<PFE> mini = ring.getPolynomial(base.getElement(2), base.getElement(0), base.getElement(1));
		FiniteField f25 = FiniteField.getFiniteField(mini, base);
		// MatrixAlgebra<PrimeFieldElement> algebra = new MatrixAlgebra<>(base, 2);
		System.exit(0);
		new TestEllipticCurve<>(f25, f25.one(), f25.one());
		new TestEllipticCurve<>(f25, f25.negative(f25.one()), f25.zero());

		mini = ring.getPolynomial(base.getElement(1), base.getElement(1), base.getElement(0), base.getElement(1));
		FiniteField f125 = FiniteField.getFiniteField(mini, base);
		EllipticCurve<FFE> ec = new EllipticCurve<>(f125, f125.one(), f125.one());

		System.out.println("Counting points!");
		System.out.println("Number of points: " + ec.getNumberOfElements());

		int n = 2;
		List<ProjectivePoint<FFE>> torsion = ec.getTorsionPoints(n);
		System.out.println(torsion.size());
		for (ProjectivePoint<FFE> point1 : torsion) {
			for (ProjectivePoint<FFE> point2 : torsion) {
				System.out.println("P1: " + point1);
				System.out.println("P2: " + point2);
				BigInteger nBig = BigInteger.valueOf(n);
				System.out.println("Weil(P1, P2): " + ec.weilPairing(nBig, point1, point2));
				System.out.println("Weil(2P1, P2): " + ec.weilPairing(nBig, ec.multiply(2, point1), point2));
				System.out.println("Weil(P1, 2P2): " + ec.weilPairing(nBig, point1, ec.multiply(2, point2)));
			}
		}
		for (int i = 0; i < 20; i++) {
			FFE re = f125.getRandomElement();
			if (!re.equals(f125.zero())) {
//				CoordinateRingElement<PFE> rp = coord.getEmbedding(f125.asGenericExtensionFieldOverPrime()
//						.asPolynomial(f125.asGenericPrimeExtensionFieldElement(re)));
//				System.out.println(re);
//				System.out.println(rp);
//				System.out.println(coord.isUnit(rp));
//				System.out.println(f125.inverse(re));
//				System.out.println(coord.inverse(rp));

				ProjectivePoint<FFE> point1 = ec.getRandomElement();
				BigInteger order1 = ec.getOrder(point1);
				ProjectivePoint<FFE> point2 = ec.getRandomElement();
				BigInteger order2 = ec.getOrder(point2);
				ProjectivePoint<FFE> point3 = ec.getRandomElement();
				BigInteger order3 = ec.getOrder(point3);
				ProjectivePoint<FFE> point4 = ec.getRandomElement();
				BigInteger order4 = ec.getOrder(point4);
				ProjectivePoint<FFE> point5 = ec.negative(ec.add(point2, point3));
				BigInteger order5 = ec.getOrder(point5);
				System.out.println("P1: " + point1 + " order: " + order1);
				System.out.println("P2: " + point2 + " order: " + order2);
				System.out.println("P3: " + point3 + " order: " + order3);
				System.out.println("P4: " + point4 + " order: " + order4);
				System.out.println("-(P2+P3): " + point5 + " order: " + order5);
				RationalFunction<FFE> f = ec.getRationalFunction(point1, point2, point3);
				System.out.println(f);
//				f.simplify();
//				System.out.println(f);
				FunctionField<FFE> ff = ec.getFunctionField();
				RationalFunction<FFE> f1 = ff.add(f, ff.one());
				RationalFunction<FFE> fsq = ff.multiply(f, f);
				System.out.println("f/f = " + ff.divide(f, f));
				System.out.println("f+f = " + ff.add(f, f));
				System.out.println("f*f = " + fsq);
				System.out.println("f+1 = " + f1);
				System.out.println("f(P1) = " + f.evaluate(point1));
				System.out.println("f(P2) = " + f.evaluate(point2));
				System.out.println("f(P3) = " + f.evaluate(point3));
				System.out.println("f(P4) = " + f.evaluate(point4));
				System.out.println("f(-(P2+P3)) = " + f.evaluate(point5));
				System.out.println("f(O) = " + f.evaluate(ec.neutral()));
				System.out.println("f+1(P1) = " + f1.evaluate(point1));
				System.out.println("f+1(P2) = " + f1.evaluate(point2));
				System.out.println("f+1(P3) = " + f1.evaluate(point3));
				System.out.println("f+1(P4) = " + f1.evaluate(point4));
				System.out.println("f+1(-(P2+P3)) = " + f1.evaluate(point5));
				System.out.println("f+1(O) = " + f1.evaluate(ec.neutral()));
				System.out.println("f*f(P1) = " + fsq.evaluate(point1));
				System.out.println("f*f(P2) = " + fsq.evaluate(point2));
				System.out.println("f*f(P3) = " + fsq.evaluate(point3));
				System.out.println("f*f(P4) = " + fsq.evaluate(point4));
				System.out.println("f*f(-(P2+P3)) = " + fsq.evaluate(point5));
				System.out.println("f*f(O) = " + fsq.evaluate(ec.neutral()));
			}
		}
		// new TestEllipticCurve<>(f125, f125.one(), f125.one());
		// new TestEllipticCurve<>(f125, f125.negative(f125.one()), f125.zero());
		/*
		 * TestAsymmetricCrypto test = new TestAsymmetricCrypto();
		 * test.doKeyAgreement(); test.doRelatedSupersingularCurveAttack(3);
		 * test.doRelatedSupersingularCurveAttack(3);
		 * test.doRelatedSupersingularCurveAttack(3);
		 * test.doRelatedSupersingularCurveAttack(3);
		 * test.doRelatedSupersingularCurveAttack(3);
		 * test.doRelatedSupersingularCurveAttack(3);
		 * test.doRelatedSupersingularCurveAttack(3);
		 */
		BigInteger one = BigInteger.ONE;
		BigInteger two = BigInteger.TWO;
		BigInteger three = BigInteger.valueOf(3);
		PrimeField f54 = PrimeField.getPrimeField(two.pow(5).multiply(three.pow(4)).subtract(one));
		UnivariatePolynomial<PFE> f54p = f54.getUnivariatePolynomialRing().getPolynomial(f54.one(), f54.zero(),
				f54.one());
		FiniteField f54sq = FiniteField.getFiniteField(f54p, f54);
		new TestSike<FFE>(f54sq);

		PrimeField f9754 = PrimeField.getPrimeField(two.pow(97).multiply(three.pow(54)).subtract(one));
		UnivariatePolynomial<PFE> f9754p = f9754.getUnivariatePolynomialRing().getPolynomial(f9754.one(), f9754.zero(),
				f9754.one());
		FiniteField f9754sq = FiniteField.getFiniteField(f9754p, f9754);
		new TestSike<FFE>(f9754sq);

		PrimeField fp434 = PrimeField.getPrimeField(two.pow(0xd8).multiply(three.pow(0x89)).subtract(one));
		UnivariatePolynomial<PFE> fp434p = fp434.getUnivariatePolynomialRing().getPolynomial(fp434.one(), fp434.zero(),
				fp434.one());
		FiniteField fp434sq = FiniteField.getFiniteField(fp434p, fp434);
		new TestSike<FFE>(fp434sq);

		// BigPrimeField f65537 = new BigPrimeField(BigInteger.valueOf(65537));
		// ExtensionField<BigPrimeFieldElement> f65537sq = new ExtensionField<>(2,
		// f65537);
		// new TestSike<ExtensionFieldElement<BigPrimeFieldElement>>(f65537sq);
		// BigPrimeField f65537 = new BigPrimeField(BigInteger.valueOf(65537));
		// ExtensionField<BigPrimeFieldElement> f65537sq = new ExtensionField<>(2,
		// f65537);
		// new TestSike<ExtensionFieldElement<BigPrimeFieldElement>>(f65537sq);
		System.exit(0);
		// new TestEllipticCurve<>(f65537, f65537.one(), f65537.one());
		// new TestEllipticCurve<>(f65537, f65537.negative(f65537.one()),
		// f65537.zero());
		PrimeField f65537 = PrimeField.getPrimeField(BigInteger.valueOf(65537));
		TestFindRelatedCurves<PFE> f65537test = new TestFindRelatedCurves<>(f65537, f65537.one(), f65537.getInteger(2),
				new ProjectivePoint<PFE>(f65537, f65537.getInteger(1), f65537.getInteger(2), f65537.one()));
		f65537test.doAttack();

		// BigInteger p256Prime = new
		// BigInteger("ffffffff00000001000000000000000000000000ffffffffffffffffffffffff",
		// 16);
		// BigInteger p256A = new BigInteger("-3");
		// BigInteger p256B = new BigInteger(
		// "41058363725152142129326129780047268409114441015993725554835256314039467401291");
		// BigPrimeField fp256 = new BigPrimeField(p256Prime);
		// BigInteger baseX = new BigInteger(
		// "48439561293906451759052585252797914202762949526041747995844080717082404635286");
		// BigInteger baseY = new BigInteger(
		// "36134250956749795798585127919587881956611106672985015071877198253568414405109");
		// ProjectivePoint<BigPrimeFieldElement> p256Generator = new
		// ProjectivePoint<>(fp256, fp256.getElement(baseX),
		// fp256.getElement(baseY), fp256.one());
		// new TestFindRelatedCurves(fp256, fp256.getElement(p256A),
		// fp256.getElement(p256B), p256Generator);
	}
}
