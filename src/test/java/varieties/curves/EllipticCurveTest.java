package varieties.curves;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Deque;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import org.junit.jupiter.api.Test;

//import edu.uci.ics.jung.graph.Graph;
//import edu.uci.ics.jung.graph.SparseMultigraph;
import fields.finitefields.FiniteField;
import fields.finitefields.FiniteField.FFE;
import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Rationals;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.Ring.FactorizationResult;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.numberfields.NumberField;
import fields.numberfields.NumberField.NFE;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.Monomial;
import fields.vectors.Vector;
import varieties.curves.elliptic.EllipticCurve;
import varieties.curves.elliptic.Frobenius;
import varieties.curves.elliptic.KernelIsogeny;
import varieties.projective.ProjectivePoint;

class EllipticCurveTest {
	private Map<Vector<IntE>, FFE> binomResult = new TreeMap<>();
	private FiniteField binomField = null;

	@Test
	void testFromPolynomial() {
		FiniteField field = FiniteField.getFiniteField(7, 4);
		PolynomialRing<FFE> ring = AbstractPolynomialRing.getPolynomialRing(field, 2, Monomial.GREVLEX);
		Map<Monomial, FFE> coeffs = new TreeMap<>();
		coeffs.put(ring.getMonomial(new int[] { 0, 2 }), field.getInteger(5));
		coeffs.put(ring.getMonomial(new int[] { 1, 1 }), field.getInteger(2));
		coeffs.put(ring.getMonomial(new int[] { 0, 1 }), field.getInteger(6));
		coeffs.put(ring.getMonomial(new int[] { 3, 0 }), field.getInteger(2));
		coeffs.put(ring.getMonomial(new int[] { 2, 0 }), field.getInteger(4));
		coeffs.put(ring.getMonomial(new int[] { 1, 0 }), field.getInteger(2));
		coeffs.put(ring.getMonomial(new int[] { 0, 0 }),
				field.add(field.getInteger(3), field.alpha(), field.power(field.alpha(), 7 * 7)));
		Polynomial<FFE> polynomial = ring.getPolynomial(coeffs);
		EllipticCurve.fromPolynomial(field, polynomial);
	}

	@Test
	void testFromJInvariantFiniteField() {
		FiniteField f = FiniteField.getFiniteField(5, 3);
		for (FFE element : f) {
			EllipticCurve<FFE> curve = EllipticCurve.fromJInvariant(f, element);
			assertEquals(element, curve.jInvariant());
		}
	}

	@Test
	void testFromJInvariantNumberField() {
		Rationals q = Rationals.q();
		NumberField f = NumberField.getNumberField(q.getUnivariatePolynomialRing().getPolynomial(q.one(), q.zero(), q.one()));
		for (int i = 0; i < 500; i++) {
			NFE element = f.getRandomElement();
			EllipticCurve<NFE> curve = EllipticCurve.fromJInvariant(f, element);
			assertEquals(element, curve.jInvariant());
		}
	}

	private Map<Integer, Polynomial<IntE>> divisionPolynomials = null;

	private <T extends Element<T>> Polynomial<T> divisionPolynomial(Field<T> field, int m) {
		Integers z = Integers.z();
		if (this.divisionPolynomials == null) {
			this.divisionPolynomials = new TreeMap<>();
		}
		if (this.divisionPolynomials.containsKey(m)) {
			PolynomialRing<T> affineRing = AbstractPolynomialRing.getPolynomialRing(field, 4, Monomial.GREVLEX);
			Polynomial<T> p = affineRing.getEmbedding(this.divisionPolynomials.get(m), new MathMap<>() {

				@Override
				public T evaluate(IntE t) {
					return field.getInteger(t);
				}
			});
			if (m % 2 == 0) {
				return affineRing.multiply(2, affineRing.getVar(2), p);
			} else {
				return p;
			}
		}
		PolynomialRing<IntE> r = AbstractPolynomialRing.getPolynomialRing(z, 4, Monomial.GREVLEX);
		Polynomial<IntE> o = r.one();
		Polynomial<IntE> x = r.getVar(1);
		Polynomial<IntE> x2 = r.getVarPower(1, 2);
		Polynomial<IntE> x3 = r.getVarPower(1, 3);
		Polynomial<IntE> x4 = r.getVarPower(1, 4);
		Polynomial<IntE> x6 = r.getVarPower(1, 6);
		Polynomial<IntE> a = r.getVar(3);
		Polynomial<IntE> a2 = r.getVarPower(3, 2);
		Polynomial<IntE> a3 = r.getVarPower(3, 3);
		Polynomial<IntE> b = r.getVar(4);
		Polynomial<IntE> b2 = r.getVarPower(4, 2);
		Polynomial<IntE> ab = r.multiply(a, b);
		Polynomial<IntE> y2 = r.add(x3, r.multiply(a, x), b);
		if (m == 0) {
			this.divisionPolynomials.put(0, r.zero());
		} else if (m == 1) {
			this.divisionPolynomials.put(1, o);
		} else if (m == 2) {
			this.divisionPolynomials.put(2, o);
		} else if (m == 3) {
			Polynomial<IntE> p = r.multiply(3, x4);
			p = r.add(p, r.multiply(6, a, x2));
			p = r.add(p, r.multiply(12, b, x));
			p = r.add(p, r.multiply(-1, a2));
			this.divisionPolynomials.put(3, p);
		} else if (m == 4) {
			Polynomial<IntE> p = x6;
			p = r.add(p, r.multiply(5, a, x4));
			p = r.add(p, r.multiply(20, b, x3));
			p = r.add(p, r.multiply(-5, a2, x2));
			p = r.add(p, r.multiply(-4, ab, x));
			p = r.add(p, r.multiply(-8, b2));
			p = r.add(p, r.multiply(-1, a3));
			p = r.multiply(2, p);
			this.divisionPolynomials.put(4, p);
		} else if (m % 2 == 0) {
			int n = m / 2;
			for (int i = n - 2; i <= n + 2; i++) {
				this.divisionPolynomial(field, i);
			}
			Polynomial<IntE> psiNm2 = this.divisionPolynomials.get(n - 2);
			Polynomial<IntE> psiNm1 = this.divisionPolynomials.get(n - 1);
			Polynomial<IntE> psiN = this.divisionPolynomials.get(n);
			Polynomial<IntE> psiN1 = this.divisionPolynomials.get(n + 1);
			Polynomial<IntE> psiN2 = this.divisionPolynomials.get(n + 2);
			Polynomial<IntE> p = r.subtract(r.multiply(psiN2, r.power(psiNm1, 2)),
					r.multiply(psiNm2, r.power(psiN1, 2)));
			p = r.multiply(psiN, p);
			this.divisionPolynomials.put(m, p);
		} else if (m % 2 == 1) {
			int n = (m - 1) / 2;
			for (int i = n - 1; i <= n + 2; i++) {
				this.divisionPolynomial(field, i);
			}
			Polynomial<IntE> psiNm1 = this.divisionPolynomials.get(n - 1);
			Polynomial<IntE> psiN = this.divisionPolynomials.get(n);
			Polynomial<IntE> psiN1 = this.divisionPolynomials.get(n + 1);
			Polynomial<IntE> psiN2 = this.divisionPolynomials.get(n + 2);
			Polynomial<IntE> firstTerm = r.multiply(psiN2, r.power(psiN, 3));
			Polynomial<IntE> secondTerm = r.multiply(psiNm1, r.power(psiN1, 3));
			if (n % 2 == 0) {
				firstTerm = r.multiply(16, firstTerm, y2, y2);
			} else {
				secondTerm = r.multiply(16, secondTerm, y2, y2);
			}
			this.divisionPolynomials.put(m, r.subtract(firstTerm, secondTerm));
		}
		return this.divisionPolynomial(field, m);

	}

	@Test
	void testDivisionPolynomial() {
		fail();
		PrimeField field = PrimeField.getPrimeField(19);
		System.out.println(divisionPolynomial(field, 19));
		Polynomial<PFE> t = divisionPolynomial(field, 19);
		PolynomialRing<PFE> polynomials = t.getPolynomialRing();
		List<PFE> list = new ArrayList<>();
		list.add(field.zero());
		list.add(field.zero());
		list.add(null);
		list.add(null);
		System.out.println(polynomials.partiallyEvaluate(t, list));
	}

	@Test
	void testSupersingularIsogenies() {
		Integers z = Integers.z();
		IntE prime = z.subtract(z.multiply(z.power(z.getInteger(2), 5), z.power(z.getInteger(3), 4)), z.one());
		FiniteField f = FiniteField.getFiniteField(prime.intValueExact(), 2);
		System.out.println(f);
		System.out.println(EllipticCurve.numberOfSupersingularCurves(prime) + " supersingular curves");
		EllipticCurve<FFE> start = EllipticCurve.fromJInvariant(f, f.getInteger(1728));
		System.out.println(start);
		System.out.println("isogenies of degree 2");
		Set<FFE> jInvariantsVisited = new TreeSet<>();
		Deque<EllipticCurve<FFE>> queue = new LinkedList<>();
		// Graph<EllipticCurve<FFE>, ProjectivePoint<FFE>> twoIsogenies = new
		// SparseMultigraph<>();
		queue.add(start);
		while (!queue.isEmpty()) {
			EllipticCurve<FFE> curve = queue.poll();
			if (jInvariantsVisited.contains(curve.jInvariant())) {
				continue;
			}
			jInvariantsVisited.add(curve.jInvariant());
			int jInvariantPower = f.smallestSubfieldContaining(curve.jInvariant());
			if (f.smallestSubfieldContaining(curve.getA()) > jInvariantPower
					|| f.smallestSubfieldContaining(curve.getB()) > jInvariantPower) {
				System.out.print("Redefining " + curve);
				curve = EllipticCurve.fromJInvariant(f, curve.jInvariant());
				System.out.println(" to " + curve);
			}
			assertTrue(curve.isSupersingular());

			List<ProjectivePoint<FFE>> twoTorsionPoints = curve.getTorsionPoints(2);
			Frobenius<FFE> frobenius = new Frobenius<>(curve, 1);
			if (!frobenius.getRange().equals(curve)) {
				// System.out.println("Skipping curve not over Fp");
			} else {
				System.out.println(curve.jInvariant());
				List<ProjectivePoint<FFE>> basis = computeBasis(2, f, curve);
				ProjectivePoint<FFE> basis1 = basis.get(0);
				ProjectivePoint<FFE> basis2 = basis.get(1);
				while (!basis1.equals(curve.neutral())) {
					System.out.println(basis1);
					System.out.println(Math.max(f.smallestSubfieldContaining(basis1.getCoord(1)),
							f.smallestSubfieldContaining(basis1.getCoord(2))));
					System.out.println(basis2);
					System.out.println(Math.max(f.smallestSubfieldContaining(basis2.getCoord(1)),
							f.smallestSubfieldContaining(basis2.getCoord(2))));
					basis1 = curve.multiply(2, basis1);
					basis2 = curve.multiply(2, basis2);
				}
			}
			List<FFE> neighbors = new ArrayList<>();
			for (ProjectivePoint<FFE> point : twoTorsionPoints) {
				if (point.equals(curve.neutral())) {
					continue;
				}
				KernelIsogeny<FFE> isogeny = new KernelIsogeny<>(curve, point, BigInteger.TWO);
				neighbors.add(isogeny.getRange().jInvariant());
				queue.add(isogeny.getRange());
				// System.out.println(curve.jInvariant() + " -> " +
				// isogeny.getRange().jInvariant());
			}
			// twoIsogenies.put(curve.jInvariant(), neighbors);
			// System.out.println();
		}
	}

	@Test
	void testSupersingularCurves() {
		fail();
		Integers z = Integers.z();
		Iterator<IntE> primes = z.primes();
		primes.next();
		primes.next();
		for (int i = 0; i < 10; i++) {
			IntE prime = primes.next();
			FiniteField f = FiniteField.getFiniteField(prime.intValueExact(), 2);
			System.out.println(f);
			System.out.println(EllipticCurve.numberOfSupersingularCurves(prime) + " supersingular curves");
			UnivariatePolynomial<FFE> hasse = hassePolynomial(f);
			System.out.println(hasse);
			System.out.println(f.factorization(hasse));
			Map<FFE, Integer> roots = f.roots(hasse);
			Set<FFE> jInvariants = new TreeSet<>();
			for (FFE root : roots.keySet()) {
				EllipticCurve<FFE> curve = EllipticCurve.fromLegendre(f, root);
				if (jInvariants.contains(curve.jInvariant())) {
					continue;
				}
				System.out.println("Curve: " + curve + " (j = " + curve.jInvariant() + ")");
				jInvariants.add(curve.jInvariant());
				assertTrue(curve.isSupersingular());
			}
		}
	}

//	@Test
//	void test() {
//		Integers z = Integers.z();
//		Iterator<IntE> it = z.primes();
//		it.next();
//		it.next();
//		for (int i = 0; i < 250; i++) {
//			binomResult.clear();
//			IntE prime = it.next();
//			FiniteField f = FiniteField.getFiniteField(prime.intValueExact(), 2);
//			System.out.println(f);
//			UnivariatePolynomial<FFE> hasse = hassePolynomial(f);
//			System.out.println(hasse);
//			System.out.println(f.factorization(hasse));
//			Map<FFE, Integer> roots = f.roots(hasse);
//			for (FFE root : roots.keySet()) {
//				EllipticCurve<FFE> curve = EllipticCurve.fromLengendre(f, root);
//				System.out.println("Curve: " + curve + " (j = " + curve.jInvariant() + ")");
//				int numberOfPoints = 0;
//				for (ProjectivePoint<FFE> point : curve) {
//					numberOfPoints++;
//				//	System.out.println(numberOfPoints + ": " + point);
//				}
//				Polynomial<FFE> divPoly = curve.getDivisionPolynomial(prime.intValueExact());
//				System.out.println(divPoly);
//				UnivariatePolynomial<FFE> divPolyRoot = f.getUnivariatePolynomialRing().getEmbedding(divPoly.getPolynomialRing().characteristicRoot(divPoly), new int[] {0});
//				System.out.println(f.factorization(divPolyRoot));
//				assertTrue(divPoly.getPolynomialRing().hasCharacteristicRoot(divPoly));
//				System.out.println(divPoly.getPolynomialRing().characteristicRoot(divPoly));
//				System.out.println("#E(" + f + ") = " + numberOfPoints + " (" + curve.getNumberOfElements() + ")");
//				assertTrue(curve.isSupersingular());
//			}
//		}
//	}

	private UnivariatePolynomial<FFE> hassePolynomial(FiniteField f) {
		UnivariatePolynomialRing<FFE> r = f.getUnivariatePolynomialRing();
		int m = (f.characteristic().intValueExact() - 1) / 2;
		List<FFE> c = new ArrayList<>();
		for (int i = 0; i <= m; i++) {
			FFE b = binom(f, m, i);
			c.add(f.multiply(b, b));
		}
		return r.getPolynomial(c);
	}

	private FFE binom(FiniteField f, int n, int k) {
		if (f != binomField) {
			binomField = f;
			binomResult.clear();
		}
		Integers z = Integers.z();
		if (k == 0) {
			return f.one();
		}
		if (k < 0 || k > n) {
			return f.zero();
		}
		IntE nInt = z.getInteger(n);
		IntE kInt = z.getInteger(k);
		Vector<IntE> pair = new Vector<>(nInt, kInt);
		if (binomResult.containsKey(pair)) {
			return binomResult.get(pair);
		}
		FFE result = f.add(binom(f, n - 1, k - 1), binom(f, n - 1, k));
		binomResult.put(pair, result);
		return result;
	}

	private List<ProjectivePoint<FFE>> computeBasis(int prime, FiniteField field, EllipticCurve<FFE> curve) {
		Integers z = Integers.z();
		IntE order = z.getInteger(curve.getNumberOfElements());
		order = z.archimedeanValue(z.sqrt(order).keySet().iterator().next());
		FactorizationResult<IntE, IntE> factors = z.uniqueFactorization(order);
		IntE basisOrder = z.power(z.getInteger(prime), factors.multiplicity(z.getInteger(prime)));
		IntE complementOrder = z.divideChecked(order, basisOrder);
		IntE lessThanOrder = z.divideChecked(basisOrder, z.getInteger(prime));
		List<ProjectivePoint<FFE>> basis = new ArrayList<>();
		ProjectivePoint<FFE> basis1;
		do {
			System.out.println("Count basis1 up");

			basis1 = curve.multiply(complementOrder, curve.getRandomElement());
		} while (!curve.multiply(basisOrder, basis1).equals(curve.neutral())
				|| curve.multiply(lessThanOrder, basis1).equals(curve.neutral()));
		ProjectivePoint<FFE> basis2;
		do {
			System.out.println("Count basis2 up");
			basis2 = curve.multiply(complementOrder, curve.getRandomElement());
		} while (!curve.multiply(basisOrder, basis2).equals(curve.neutral())
				|| curve.multiply(lessThanOrder, basis2).equals(curve.neutral())
				|| curve.weilPairing(basisOrder.getValue(), basis1, basis2).equals(field.one()));
		basis.add(basis1);
		basis.add(basis2);
		return basis;
	}

}
