package varieties.curves;

import static org.junit.jupiter.api.Assertions.fail;

import java.io.IOException;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;

import org.junit.jupiter.api.Test;

import fields.finitefields.FiniteField;
import fields.finitefields.FiniteField.FFE;
import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.interfaces.Field.Extension;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import util.Pair;
import util.graph.Edge;
import util.graph.Graph;
import varieties.curves.elliptic.EllipticCurve;
import varieties.curves.elliptic.EndomorphismRing;
import varieties.curves.elliptic.Isogeny;
import varieties.projective.ProjectivePoint;

class EndomorphismRingTest {

	@Test
	void ordinaryTest() throws IOException {
		Integers z = Integers.z();
		int p = 13;
		int l = 3;
		PrimeField fp = PrimeField.getPrimeField(p);
		Set<Pair<PFE, IntE>> jInvariants = new TreeSet<>();
		for (int a = 0; a < p; a++) {
			for (int b = 0; b < p; b++) {
				if ((4 * a * a * a + 27 * b * b) % p == 0) {
					continue;
				}
				a = 2;
				b = 6;
				EllipticCurve<PFE> curve = new EllipticCurve<>(fp, fp.getInteger(a), fp.getInteger(b));
//				System.out.println(curve);
//				System.out.println("j Invariant: " + curve.jInvariant());
				if (curve.isSupersingular() || jInvariants
						.contains(new Pair<>(curve.jInvariant(), z.getInteger(curve.getNumberOfElements())))) {
					continue;
				}
				System.out.println(curve + ": " + curve.getNumberOfElements());
				System.out.println("j Invariant: " + curve.jInvariant());
				System.out.println(curve.getEndomorphismRing());
//				curve = curve.getQuadraticTwist();
//				System.out.println(curve + ": " + curve.getNumberOfElements());
//				System.out.println("j Invariant: " + curve.jInvariant());
//				Isogeny<PFE> frobenius = curve.frobeniousEndomorphism();
//				List<ProjectivePoint<PFE>> points = curve.getTorsionPoints(3);
//				for (ProjectivePoint<PFE> point : points) {
//					System.out.println(curve.add(
//							curve.add(frobenius.evaluate(frobenius.evaluate(point)),
//									curve.multiply(curve.trace().negate(), frobenius.evaluate(point))),
//							curve.multiply(p % 3, point)));
//				}
//				System.out.println();
//				Extension<PFE, PFE, FFE, FiniteField> ext = fp
//						.getExtension(FiniteField.getFiniteField(fp, 6).minimalPolynomial());
//				EllipticCurve<FFE> extCurve = curve.extendBaseField(ext).getCurve();
//				Isogeny<FFE> extFrobenius = extCurve.frobeniousEndomorphism();
//				List<ProjectivePoint<FFE>> extPoints = extCurve.getTorsionPoints(3);
//				// 3*X^3 + 6*X^2 + -4*X + -2
//				UnivariatePolynomialRing<FFE> ring = ext.extension().getUnivariatePolynomialRing();
//				UnivariatePolynomial<FFE> lhs = ring.parse("3*X^3 + -6*X^2 + -4*X + 2");
//				for (ProjectivePoint<FFE> point : extPoints) {
//					if (point.equals(extCurve.neutral())) {
//						continue;
//					}
//					FFE x = point.getDehomogenisedCoord(1, 3);
//					FFE result = ring.evaluate(lhs, x);
//					System.out.println(result);
//					System.out.println(extCurve.add(extFrobenius.evaluate(extFrobenius.evaluate(point)),
//							extCurve.multiply(ext.extension().getNumberOfElements().intValueExact() % 3, point)));
////					System.out
////							.println(
////									extCurve.add(
////											extCurve.add(extFrobenius.evaluate(extFrobenius.evaluate(point)),
////													extCurve.multiply(curve.trace().negate(),
////															extFrobenius.evaluate(point))),
////											extCurve.multiply(p % 3, point)));
//				}
//				fail();
//				Graph<EllipticCurve<PFE>, Isogeny<PFE>> isogenyGraph = curve.isogenyGraph(l);
//				for (EllipticCurve<PFE> node : isogenyGraph.vertexSet()) {
//					System.out.print("  " + node.jInvariant() + " (" + node + "): ");
//					List<String> edges = isogenyGraph.edges(node).stream()
//							.map((Edge<EllipticCurve<PFE>, Isogeny<PFE>> e) -> e.otherVertex(node).jInvariant()
//									.toString() + " (" + e.otherVertex(node).toString() + ")")
//							.collect(Collectors.toList());
//					System.out.println(String.join(", ", edges));
//					jInvariants.add(new Pair<>(node.jInvariant(), z.getInteger(node.getNumberOfElements())));
//				}
//				System.out.println();
			}
		}
	}

	@Test
	void superSingularTest() throws IOException {
//		PrimeField fp = PrimeField.getPrimeField(257);
//		PolynomialRing<PFE> polynomials = AbstractPolynomialRing.getPolynomialRing(fp, Monomial.LEX,
//				new String[] { "Y", "X", "a1", "a2", "a3", "a4", "a6", "x_r", "y_r", "s" });
//		List<Polynomial<PFE>> constraints = new ArrayList<>();
//		constraints.add(polynomials.parse("-1*Y^2 + -1*a1*X*Y + -1*a3*Y + X^3 + a2*X^2 + a4*X + a6"));
//		constraints.add(polynomials.parse("-1*y_r^2 + -1*a1*x_r*y_r + -1*a3*y_r + x_r^3 + a2*x_r^2 + a4*x_r + a6"));
//		constraints.add(polynomials.parse("Y + y_r + s*X + -1*s*x_r"));
//		PolynomialIdeal<PFE> ideal = polynomials.getIdeal(constraints);
//		for (Polynomial<PFE> polynomial : ideal.generators()) {
//			System.out.println(polynomial);
//		}
//		Polynomial<PFE> toSolve = ideal.generators().get(ideal.generators().size() - 2);
//				QuotientAndRemainderResult<Polynomial<PFE>> qr = polynomials.quotientAndRemainder(toSolve,
//				polynomials.parse("X + -1*x_r"));
//		System.out.println(qr.getQuotient());
//		System.out.println();
//		System.out.println(qr.getRemainder());
//		System.out.println();
//			UnivariatePolynomial<Polynomial<PFE>> asuni = polynomials.asUnivariatePolynomial(qr.getQuotient(), 2);
//		UnivariatePolynomialRing<Polynomial<PFE>> asuniring = asuni.getPolynomialRing();
//		Polynomial<PFE> discriminant = asuniring.discriminant(asuni);
//		UnivariatePolynomial<Polynomial<PFE>> discuni = discriminant.getPolynomialRing()
//				.asUnivariatePolynomial(discriminant, 9);
//		System.out.println(discuni);
//		fail();
//		Rationals q = Rationals.q();
//		FiniteField fq = FiniteField.getFiniteField(23*23);
//		FFE sqrt = fq.sqrt(fq.getInteger(-1)).keySet().iterator().next();
//		assertTrue(fq.multiply(sqrt, sqrt).equals(fq.getInteger(-1)));
//		EllipticCurve<FFE> curve = new EllipticCurve<>(fq, fq.getInteger(1), fq.zero());
//		Isomorphism<FFE> iota = new Isomorphism<>(curve, sqrt, fq.zero(), fq.zero(), fq.zero());
//		assertTrue(iota.getRange().equals(curve));
//		System.out.println(iota);
//		Frobenius<FFE> pi = new Frobenius<>(curve, 1);
//		assertTrue(pi.getRange().equals(curve));
//		System.out.println(pi);
//		AdditionIsogeny<FFE> added = new AdditionIsogeny<>(iota, pi);
//		assertTrue(added.getRange().equals(curve));
//		System.out.println(added);
//		List<ProjectivePoint<FFE>> kernel = new ArrayList<>();
//		for (ProjectivePoint<FFE> point : curve) {
////			System.out.println(point + ", order: " + curve.getOrder(point));
////			System.out.println(iota.evaluate(point) + " + " + pi.evaluate(point) + " = " + added.evaluate(point));
////			System.out.println();
//			if (added.evaluate(point).equals(curve.neutral())) {
//				kernel.add(point);
//			}
//		}
//		// System.out.println(kernel);
//		// fail();
//		Set<ProjectivePoint<FFE>> halfKernel = new TreeSet<>();
//		for (ProjectivePoint<FFE> kernelPoint : kernel) {
//			halfKernel.add(curve.multiply(2, kernelPoint));
//		}
//		// System.out.println(halfKernel);
//		KernelIsogeny<FFE> halved = new KernelIsogeny<>(curve, halfKernel);
//		// System.out.println(halved.getRange());
//		List<Isogeny<FFE>> isogenies = new ArrayList<>();
//		for (Isomorphism<FFE> isomorphism : halved.getRange().getIsomorphisms(curve)) {
//			isogenies.add(new CompositionIsogeny<FiniteField.FFE>(halved, isomorphism));
//		}
//		System.out.println();
//		Isogeny<FFE> result = null;
//		for (Isogeny<FFE> isogeny : isogenies) {
//			boolean success = true;
//			for (ProjectivePoint<FFE> point : curve) {
//				// System.out.println(point + ", order: " + curve.getOrder(point));
//				ProjectivePoint<FFE> eval = isogeny.evaluate(point);
//				// System.out.println("2*" + eval + " = " + curve.multiply(2, eval) + " = " +
//				// added.evaluate(point));
//				if (!curve.multiply(2, eval).equals(added.evaluate(point))) {
//					success = false;
//					break;
//				}
////				System.out.println();
//			}
//			if (success) {
//				result = isogeny;
//			}
//		}
//		System.out.println(result);
//		System.out.println(curve.getNumberOfElements());
//		System.out.println(halfKernel.size());
//		List<ProjectivePoint<FFE>> addedEval = new ArrayList<>();
//		for (ProjectivePoint<FFE> basis : curve.getTorsionPointBasis((fq.characteristic().intValueExact() + 1) / 2)) {
//			addedEval.add(added.evaluate(basis));
//		}
//		System.out.println(addedEval);
//		System.out.println(curve.weilPairing(addedEval.get(0), addedEval.get(1)));
//		Set<ProjectivePoint<FFE>> asSet = new TreeSet<>();
//		asSet.addAll(addedEval);
//		System.out.println(new KernelIsogeny<>(curve, asSet));
//		fail();
//		RationalQuaternions quaternions = RationalQuaternions.quaternions(q.getInteger(-1), q.getInteger(-23));
//		RationalQuaternionOrder order = new RationalQuaternionOrder(quaternions, quaternions.i(),
//				quaternions.divide(quaternions.add(quaternions.i(), quaternions.j()), quaternions.getInteger(2)),
//				quaternions.divide(quaternions.add(quaternions.one(), quaternions.k()), quaternions.getInteger(2)));
//		System.out.println(order.isMaximal());
//		System.out.println(quaternions.reducedNorm(order.getModuleGenerators().get(1)));
//		System.out.println(quaternions.reducedNorm(order.getModuleGenerators().get(2)));
//		System.out.println(quaternions.reducedNorm(order.getModuleGenerators().get(3)));
//		System.out.println(quaternions.power(order.getModuleGenerators().get(1), 2));
//		System.out.println(quaternions.power(order.getModuleGenerators().get(2), 2));
//		System.out.println(quaternions.power(order.getModuleGenerators().get(3), 2));
		// p = 3 mod 4
		// i^2 = -1
		// j^2 = -p
		// ij = -ji
		// ij = k
		// k^2 = -p
		// i*k = i*i*j = -j
		// k*i = i*j*i = j
		// j*k = j*i*j = p*i
		// k*j = i*j*j = -p*i
		// 1, i, (i + j)/2, (1+k)/2
		// ((i+j)/2)^2 = 1/4*(i^2 + j^2) = -1/4*(p+1)
		FiniteField field = FiniteField.getFiniteField(43 * 43);
		EllipticCurve<FFE> superSingular = new EllipticCurve<>(field, field.getInteger(-1), field.zero());
		EllipticCurve<FFE> ordinary = EllipticCurve.fromJInvariant(field, field.zero());
		EndomorphismRing<FFE> superSingularEndomorphisms = superSingular.getEndomorphismRing();
		EndomorphismRing<FFE> ordinaryEndomorphisms = ordinary.getEndomorphismRing();
	}

}
