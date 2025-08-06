package cryptography;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.nio.charset.Charset;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.junit.jupiter.api.Test;

import cryptography.UnbalancedOilAndVinegar.UOVPrivateKey;
import fields.finitefields.Char2Field;
import fields.finitefields.Char2Field.C2E;
import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.helper.TranscendentalFieldExtension;
import fields.helper.TranscendentalFieldExtension.TExt;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.CoordinateRing;
import fields.polynomials.Monomial;
import fields.polynomials.PolynomialIdeal;
import fields.vectors.Matrix;
import fields.vectors.MatrixModule;
import fields.vectors.Vector;
import util.Pair;
import varieties.affine.AffinePoint;
import varieties.affine.AffineScheme;
import varieties.affine.AffineScheme.BlowUpResult;

class UnbalancedOilAndVinegarTest {

	@Test
	void test() {
		UnbalancedOilAndVinegar uov = new UnbalancedOilAndVinegar(4, 10, 3, Sha3.SHAKE_256);
		byte[] msg = "This is a test message".getBytes(Charset.forName("UTF8"));
		UOVPrivateKey privateKey = uov.createPrivateKey();
		Vector<Polynomial<C2E>> publicKey = uov.createPublicKey(privateKey);
		Vector<C2E> signature = uov.sign(msg, privateKey);
		assertTrue(uov.verify(msg, signature, publicKey));
	}

	@Test
	void forgeryTest() {
		UnbalancedOilAndVinegar uov = new UnbalancedOilAndVinegar(2, 4, 3, Sha3.SHAKE_256);
		byte[] msg = "This is a test message".getBytes(Charset.forName("UTF8"));
		UOVPrivateKey privateKey = uov.createPrivateKey();
		Vector<Polynomial<C2E>> publicKey = uov.createPublicKey(privateKey);
		Vector<C2E> signature = uov.sign(msg, privateKey);
		assertTrue(uov.verify(msg, signature, publicKey));
		byte[] secondMsg = "This is a second test message".getBytes(Charset.forName("UTF8"));
		Vector<C2E> parsed = uov.parse(secondMsg);
		Vector<C2E> forgedSignature = uov.getPolynomialRing().solve(publicKey, parsed, 1).get(0);
		assertTrue(uov.verify(secondMsg, forgedSignature, publicKey));
	}

	<T extends Element<T>> Matrix<T> toMatrix(Field<T> field, PolynomialRing<T> polynomials, Polynomial<T> t) {
		if (t.degree() != 2 || !polynomials.isHomogeneous(t)) {
			throw new ArithmeticException("not a quadratic form");
		}
		List<List<T>> result = new ArrayList<>();
		for (int i = 0; i < polynomials.numberOfVariables(); i++) {
			List<T> row = new ArrayList<>();
			for (int j = 0; j < polynomials.numberOfVariables(); j++) {
				if (j < i) {
					row.add(field.zero());
				} else {
					int[] exp = new int[polynomials.numberOfVariables()];
					exp[i]++;
					exp[j]++;
					row.add(t.coefficient(polynomials.getMonomial(exp)));
				}
			}
			result.add(row);
		}
		return new Matrix<>(result);
	}

	<T extends Element<T>> Polynomial<T> toPolynomial(Field<T> field, PolynomialRing<T> polynomials, Matrix<T> t) {
		Map<Monomial, T> result = new TreeMap<>();
		for (int i = 0; i < polynomials.numberOfVariables(); i++) {
			for (int j = i; j < polynomials.numberOfVariables(); j++) {
				T coeff = t.entry(i + 1, j + 1);
				if (i != j) {
					coeff = field.add(t.entry(j + 1, i + 1), coeff);
				}
				int[] exp = new int[polynomials.numberOfVariables()];
				exp[i]++;
				exp[j]++;
				result.put(polynomials.getMonomial(exp), coeff);
			}
		}
		return polynomials.getPolynomial(result);
	}

	// @Test
	void blowupTest() {
		int oil = 2;
		int vinegar = 3;
		UnbalancedOilAndVinegar uov = new UnbalancedOilAndVinegar(oil, vinegar, 4, Sha3.SHAKE_256);
		UOVPrivateKey privateKey = uov.createPrivateKey();
		Vector<Polynomial<C2E>> publicKey = uov.createPublicKey(privateKey);
		AffineScheme<C2E> uovScheme = new AffineScheme<>(uov.getField(),
				uov.getPolynomialRing().getIdeal(publicKey.asList()).divideOut());
		List<C2E> originList = new ArrayList<>();
		for (int i = 0; i < oil + vinegar; i++) {
			originList.add(uov.getField().zero());
		}
		AffinePoint<C2E> origin = new AffinePoint<>(uov.getField(), originList);
		BlowUpResult<C2E> blowup = uovScheme.blowupPoint(origin);
		System.out.println(blowup.getBlowUp().getDomain().getAffineCover().get(0));
		for (int i = 0; i < blowup.getBlowUp().getDomain().getAffineCover().size(); i++) {
			System.out.println(blowup.getBlowUp().getDomain().getAffineCover().get(i).getCotangentSyzygies());
			System.out.println(
					blowup.getBlowUp().getDomain().getAffineCover().get(i).singularLocus().get().simplifyDomain());
		}
		System.out.println(blowup.getInverse().getMorphism().closure().simplifyDomain());
	}
	/*
	 * @Test void flatnessTest() { int oil = 2; int vinegar = 3;
	 * UnbalancedOilAndVinegar uov = new UnbalancedOilAndVinegar(oil, vinegar, 4,
	 * Sha3.SHAKE_256); UOVPrivateKey privateKey = uov.createPrivateKey();
	 * AffineScheme<C2E> domain = AffineScheme.affineSpace(uov.getField(), oil +
	 * vinegar); AffineScheme<C2E> range = AffineScheme.affineSpace(uov.getField(),
	 * oil); AffineMorphism<C2E> privateMorphism =
	 * AffineMorphism.fromPolynomials(domain, range,
	 * privateKey.getPolynomials().asList()); boolean found = false;
	 * MatrixModule<C2E> mm = new MatrixModule<>(uov.getField(), oil, oil);
	 * AffinePoint<C2E> point; do { List<C2E> eval = new ArrayList<>(); for (int i =
	 * 0; i < vinegar; i++) { eval.add(uov.getField().getRandomElement()); } for
	 * (int i = 0; i < oil; i++) { eval.add(null); } List<List<C2E>> rows = new
	 * ArrayList<>(); for (int i = 0; i < oil; i++) { Polynomial<C2E> evaluated =
	 * uov.getPolynomialRing() .partiallyEvaluate(privateKey.getPolynomials().get(i
	 * + 1), eval); List<C2E> row = new ArrayList<>(); for (int j = 0; j < oil; j++)
	 * { int[] exp = new int[oil + vinegar]; exp[vinegar + j] = 1;
	 * row.add(evaluated.coefficient(uov.getPolynomialRing().getMonomial(exp))); }
	 * rows.add(row); } List<C2E> pointList = new ArrayList<>();
	 * pointList.addAll(eval.subList(0, vinegar)); for (int i = 0; i < oil; i++) {
	 * pointList.add(uov.getField().getRandomElement()); } point = new
	 * AffinePoint<>(uov.getField(), pointList); Matrix<C2E> matrix = new
	 * Matrix<>(rows); found = mm.rank(matrix) != oil; if (found) {
	 * System.out.println(point); List<Vector<C2E>> kernel = mm.kernelBasis(matrix);
	 * System.out.println(kernel);
	 * System.out.println(privateMorphism.evaluate(point)); List<C2E> tangent = new
	 * ArrayList<>(); for (int i = 0; i < vinegar; i++) {
	 * tangent.add(point.getCoord(i + 1)); } for (int i = 0; i < oil; i++) {
	 * tangent.add(uov.getField().add(point.getCoord(vinegar + i + 1),
	 * kernel.get(0).get(i + 1))); } System.out.println(tangent);
	 * System.out.println(privateMorphism.evaluate(new AffinePoint<>(uov.getField(),
	 * tangent))); } } while (!found); AffinePoint<C2E> evaluated =
	 * privateMorphism.evaluate(point); DifferentialForms<C2E> omega =
	 * uov.getPolynomialRing().differentialForms(); ExteriorProduct<C2E,
	 * Vector<C2E>> product = omega.evaluationSpace(); List<Vector<C2E>>
	 * cotangentBasis = new ArrayList<>(); for (int i = 0; i < oil; i++) {
	 * cotangentBasis.add(product.asGradedVector(
	 * omega.evaluate(omega.derivative(privateKey.getPolynomials().get(i + 1)),
	 * point), 1)); } PolynomialRing<C2E> eliminationRing =
	 * uov.getPolynomialRing().addVariableWithElimination(2); List<Polynomial<C2E>>
	 * eliminationBasis = new ArrayList<>(); for (Vector<C2E> cotangent :
	 * cotangentBasis) { eliminationBasis.add(eliminationRing
	 * .getEmbeddingWithElimination(uov.getPolynomialRing().getLinear(cotangent.
	 * asList()), 2)); } eliminationBasis.add(eliminationRing.subtract(
	 * eliminationRing.multiply(eliminationRing.getVar(1),
	 * eliminationRing.getVar(2)), eliminationRing.one())); List<Polynomial<C2E>>
	 * subs = new ArrayList<>(); for (int i = 0; i < oil + vinegar; i++) {
	 * subs.add(eliminationRing.add(
	 * eliminationRing.getEmbeddingWithElimination(uov.getPolynomialRing().getVar(i
	 * + 1), 2), eliminationRing.multiply(point.getCoord(i + 1),
	 * eliminationRing.getVar(1)))); } System.out.println(point);
	 * System.out.println(subs); System.out.println(evaluated); for (int i = 0; i <
	 * oil; i++) { System.out.println(privateKey.getPolynomials().get(i + 1));
	 * eliminationBasis.add(
	 * eliminationRing.subtract(eliminationRing.substitute(privateKey.getPolynomials
	 * ().get(i + 1), subs), eliminationRing.multiply(evaluated.getCoord(i + 1),
	 * eliminationRing.getVarPower(1, 2)))); } for (Polynomial<C2E> generator :
	 * eliminationBasis) { System.out.println(generator); } PolynomialIdeal<C2E>
	 * ideal = eliminationRing.getIdeal(eliminationBasis);
	 * System.out.println(ideal.dimension()); for (Polynomial<C2E> generator :
	 * ideal.generators()) { // if (generator.degree(1) <= 0) {
	 * System.out.println(generator); // } }
	 * System.out.println(eliminationRing.solve(ideal).size());
	 * System.out.println(); AffinePoint<C2E> maxPoint = null; int max = -1;
	 * AffinePoint<C2E> minPoint = null; int min = -1; for (int i = 0; i < 20; i++)
	 * { AffinePoint<C2E> rand = domain.getRandomElement(); System.out.println();
	 * for (int j = 0; j < oil; j++) { System.out.println(product.asGradedVector(
	 * omega.evaluate(omega.derivative(privateKey.getPolynomials().get(j + 1)),
	 * rand), 1)); } // AffineMorphism<C2E> preimage = privateMorphism //
	 * .preimage(range.pointAsMorphism(rand).restrict(SpectrumOfField.POINT).
	 * getMorphism()); // System.out.println(preimage.getDomain().dimension() + " "
	 * + preimage.getDomain().degree()); // List<AffinePoint<C2E>> list =
	 * privateMorphism.preimageList(rand); // System.out.println(list.size()); // if
	 * (min == -1 || list.size() < min) { // minPoint = rand; // min = list.size();
	 * // } // if (max == -1 || list.size() > min) { // maxPoint = rand; // max =
	 * list.size(); // } } System.out.println(minPoint);
	 * System.out.println(maxPoint); }
	 */
	/*
	 * @Test void relativeDifferentialFormsTest() { int oil = 4; int vinegar = 7;
	 * UnbalancedOilAndVinegar uov = new UnbalancedOilAndVinegar(oil, vinegar, 8,
	 * Sha3.SHAKE_256); UOVPrivateKey privateKey = uov.createPrivateKey();
	 * Vector<Polynomial<C2E>> publicKey = uov.createPublicKey(privateKey); byte[]
	 * message = "This is a test message".getBytes(Charset.forName("UTF8"));
	 * Vector<C2E> signature = uov.sign(message, privateKey);
	 * assertTrue(uov.verify(message, signature, publicKey)); AffineScheme<C2E>
	 * domain = AffineScheme.affineSpace(uov.getField(), oil + vinegar);
	 * AffineScheme<C2E> range = AffineScheme.affineSpace(uov.getField(), oil);
	 * AffineMorphism<C2E> privateMorphism = AffineMorphism.fromPolynomials(domain,
	 * range, privateKey.getPolynomials().asList()); RelativeDifferentialForms<C2E>
	 * relativeDifferentialForms = privateMorphism.getRelativeDifferentialForms();
	 * for (int i = 1; i <= vinegar; i++) {
	 * System.out.println(relativeDifferentialForms.getSyzygies(i)); } }
	 */

	@Test
	void universalSignerTest() {
		int oil = 4;
		int vinegar = 7;
		UnbalancedOilAndVinegar uov = new UnbalancedOilAndVinegar(oil, vinegar, 8, Sha3.SHAKE_256);
		UOVPrivateKey privateKey = uov.createPrivateKey();
		Vector<Polynomial<C2E>> publicKey = uov.createPublicKey(privateKey);
		Char2Field field = uov.getField();
		TranscendentalFieldExtension<C2E> base = new TranscendentalFieldExtension<>(field, vinegar);
		PolynomialRing<TExt<C2E>> functionField = AbstractPolynomialRing.getPolynomialRing(base, oil, Monomial.GREVLEX);
		PolynomialRing<C2E> polynomialRing = AbstractPolynomialRing.getPolynomialRing(field, oil + vinegar,
				Monomial.GREVLEX);
		List<Polynomial<C2E>> baseChange = new ArrayList<>();
		for (int i = 0; i < vinegar; i++) {
			Polynomial<C2E> changed = polynomialRing.getVar(i + 1);
			for (int j = 0; j < oil; j++) {
				changed = polynomialRing.add(polynomialRing.multiply(privateKey.getOilMatrix().entry(i + 1, j + 1),
						polynomialRing.getVar(vinegar + j + 1)), changed);
			}
			baseChange.add(changed);
		}
		for (int i = 0; i < oil; i++) {
			baseChange.add(polynomialRing.getVar(vinegar + i + 1));
		}
		List<Polynomial<C2E>> privatePolynomials = new ArrayList<>();
		for (int k = 0; k < oil; k++) {
			Polynomial<C2E> polynomial = polynomialRing.substitute(publicKey.get(k + 1), baseChange);
			System.out.println(polynomial);
			privatePolynomials.add(polynomial);
		}
	}

	// @Test
	void functionFieldTest() {
		int oil = 4;
		int vinegar = 7;
//		int counter = 0;
		UnbalancedOilAndVinegar uov = new UnbalancedOilAndVinegar(oil, vinegar, 8, Sha3.SHAKE_256);
		PolynomialRing<C2E> polynomials = uov.getPolynomialRing();
		UOVPrivateKey privateKey = uov.createPrivateKey();
		Vector<Polynomial<C2E>> publicKey = uov.createPublicKey(privateKey);
		byte[] message = "This is a test message".getBytes(Charset.forName("UTF8"));
		Vector<C2E> signature = uov.sign(message, privateKey);
		assertTrue(uov.verify(message, signature, publicKey));
		Vector<C2E> parsed = uov.parse(message);
//		PolynomialIdeal<C2E> ideal1 = polynomials.getIdeal(publicKey.asList());
//		System.out.println(ideal1);
		PolynomialRing<C2E> dbl = AbstractPolynomialRing.getPolynomialRing(uov.getField(), 2 * oil + 2 * vinegar,
				Monomial.GREVLEX);
		List<Polynomial<C2E>> newPublic = new ArrayList<>();
		int[] map = new int[oil + vinegar];
		List<Polynomial<C2E>> sum = new ArrayList<>();
		for (int i = 0; i < oil + vinegar; i++) {
			map[i] = oil + vinegar + i;
			sum.add(dbl.add(dbl.getVar(i + 1), dbl.getVar(oil + vinegar + i + 1)));
		}
		for (Polynomial<C2E> polynomial : publicKey.asList()) {
			Polynomial<C2E> firstEmbedding = dbl.getEmbedding(polynomial);
			newPublic.add(firstEmbedding);
			Polynomial<C2E> secondEmbedding = dbl.getEmbedding(polynomial, map);
			newPublic.add(secondEmbedding);
			Polynomial<C2E> third = dbl.subtract(dbl.substitute(polynomial, sum),
					dbl.add(firstEmbedding, secondEmbedding));
			newPublic.add(third);
			System.out.println(firstEmbedding);
			System.out.println(secondEmbedding);
			System.out.println(third);
			int count = 0;
			for (int i = 0; i < 2 * oil + 2 * vinegar; i++) {
				for (int j = 0; j < 2 * oil + 2 * vinegar; j++) {
					int[] exp = new int[2 * oil + 2 * vinegar];
					exp[i]++;
					exp[j]++;
					if (third.coefficient(dbl.getMonomial(exp)).equals(uov.getField().zero())) {
						// System.out.println("0*X_" + i + "*X_" + j);
						count++;
					}
				}
			}
			System.out.println("uov " + count);
			// (a b) (c d) (a.c a.d)
			// (0 0).(0 0) = ( 0 0)
			// (c d) (a b) (c.a c.b)
			// (0 0).(0 0) = ( 0 0)
		}
		List<Polynomial<C2E>> random = new ArrayList<>();
		for (int k = 0; k < oil; k++) {
			Map<Monomial, C2E> poly = new TreeMap<>();
			for (int i = 0; i < oil + vinegar; i++) {
				for (int j = 0; j < oil + vinegar; j++) {
					int[] exp = new int[oil + vinegar];
					exp[i]++;
					exp[j]++;
					poly.put(polynomials.getMonomial(exp), uov.getField().getRandomElement());
				}
			}
			random.add(polynomials.getPolynomial(poly));
		}
		for (Polynomial<C2E> polynomial : random) {
			Polynomial<C2E> firstEmbedding = dbl.getEmbedding(polynomial);
			newPublic.add(firstEmbedding);
			Polynomial<C2E> secondEmbedding = dbl.getEmbedding(polynomial, map);
			newPublic.add(secondEmbedding);
			Polynomial<C2E> third = dbl.subtract(dbl.substitute(polynomial, sum),
					dbl.add(firstEmbedding, secondEmbedding));
			newPublic.add(third);
			System.out.println(firstEmbedding);
			System.out.println(secondEmbedding);
			System.out.println(third);
			int count = 0;
			for (int i = 0; i < 2 * oil + 2 * vinegar; i++) {
				for (int j = 0; j < 2 * oil + 2 * vinegar; j++) {
					int[] exp = new int[2 * oil + 2 * vinegar];
					exp[i]++;
					exp[j]++;
					if (third.coefficient(dbl.getMonomial(exp)).equals(uov.getField().zero())) {
						// System.out.println("0*X_" + i + "*X_" + j);
						count++;
					}
				}
			}
			System.out.println("random " + count);
		}
//		PolynomialIdeal<C2E> ideal2 = dbl.getIdeal(newPublic);
//		System.out.println(ideal2);
//		DifferentialForms<C2E> differentialForms = polynomials.differentialForms();
//		FiniteVectorSpace<C2E> space = new FiniteVectorSpace<>(uov.getField(), oil + vinegar);
//		List<Vector<Polynomial<C2E>>> linearEqs = new ArrayList<>();
//		for (Polynomial<C2E> generator : publicKey.asList()) {
//			DifferentialForm<C2E> differential = differentialForms.derivative(generator);
//			System.out.println(differential);
//			Vector<Polynomial<C2E>> asVector = differentialForms.asGradedVector(differential, 1);
//			linearEqs.add(asVector);
//		}
//		Matrix<Polynomial<C2E>> matrix = Matrix.fromColumns(linearEqs);
//		MatrixModule<Polynomial<C2E>> mm = matrix.getModule(polynomials);
//		PolynomialIdeal<C2E> ideal = polynomials.getIdeal(mm.minors(matrix, oil));
//		for (Polynomial<C2E> generator : ideal.generators()) {
//			System.out.println(generator);
//		}
//		System.out.println(ideal);
//		polynomials.solve(ideal, 10);
//		fail();
//	AffineScheme<C2E> domain = new AffineScheme<>(uov.getField(), AbstractPolynomialRing
//				.getPolynomialRing(uov.getField(), oil + vinegar, Monomial.GREVLEX).getZeroIdeal().divideOut());
//		AffineScheme<C2E> range = new AffineScheme<>(uov.getField(), AbstractPolynomialRing
//				.getPolynomialRing(uov.getField(), oil, Monomial.GREVLEX).getZeroIdeal().divideOut());
//		List<CoordinateRingElement<C2E>> map = new ArrayList<>();
//		for (Polynomial<C2E> polynomial : privateKey.getPolynomials().asList()) {
//			map.add(domain.getCoordinateRing().getEmbedding(polynomial));
//		}
//		AffineMorphism<C2E> morphism = new AffineMorphism<>(domain, range, map);
//		System.out.println(morphism.isSeparated());
//		Vector<C2E> parsed = uov.parse(message);
//		AffinePoint<C2E> asPoint = new AffinePoint<>(uov.getField(), parsed.asList());
//		AffineMorphism<C2E> preimage = morphism.preimage(asPoint);
//		AffinePoint<C2E> randomPoint = preimage.getDomain().getRandomElement();
//		AffinePoint<C2E> preimagePoint = preimage.evaluate(randomPoint);
//		Vector<C2E> secondSignature = new Vector<>(preimagePoint.getCoords());
//		assertTrue(uov.verify(message, secondSignature, publicKey));
//		AffineFiberedProduct<C2E> graph = morphism.graph();
//		AffineScheme<C2E> birationalDomain = graph.getProduct();
//		List<CoordinateRingElement<C2E>> birationalMap = new ArrayList<>();
//		for (int i = 0; i < vinegar; i++) {
//			birationalMap.add(birationalDomain.getCoordinateRing().getVar(i + 1));
//		}
//		for (int i = 0; i < oil; i++) {
//			birationalMap.add(birationalDomain.getCoordinateRing().getVar(oil + vinegar + i + 1));
//		}
//		AffineMorphism<C2E> birationalMorphism = new AffineMorphism<>(birationalDomain, domain, birationalMap);
//		GeneralRationalFunction<C2E, AffinePoint<C2E>, AffinePoint<C2E>> inverse = birationalMorphism.birationalInverse();
//		Optional<AffinePoint<C2E>> evaluated;
//		do {
//		List<C2E> saltAndImage = new ArrayList<>();
//		for (int i = 0; i < vinegar; i++) {
//			saltAndImage.add(uov.getField().getRandomElement());
//		}
//		saltAndImage.addAll(parsed.asList());
//		evaluated = inverse.evaluate(new AffinePoint<>(uov.getField(), saltAndImage));
//		} while (evaluated.isEmpty());
//		Vector<C2E> thirdSignature = new Vector<>(graph.getProjection1().evaluate( evaluated.get()).getCoords());
//		assertTrue(uov.verify(message, thirdSignature, publicKey));
//		System.out.println("Success");

//		TranscendentalFieldExtension<C2E> ff = new TranscendentalFieldExtension<>(uov.getField(), polynomials);
//		List<List<TExt<C2E>>> extractedLinear = new ArrayList<>();
//		List<TExt<C2E>> extractedConstant = new ArrayList<>();
//		for (int i = 0; i < oil; i++) {
//			Polynomial<C2E> polynomial = privateKey.getPolynomials().get(i + 1);
//			Map<Monomial, C2E> constant = new TreeMap<>();
//			List<Map<Monomial, C2E>> linear = new ArrayList<>();
//			for (int j = 0; j < oil; j++) {
//				linear.add(new TreeMap<>());
//			}
//			for (Monomial m : polynomial.monomials()) {
//				int oilCoeff = -1;
//				for (int j = 0; j < oil; j++) {
//					if (m.exponents()[vinegar + j] != 0) {
//						oilCoeff = j;
//						break;
//					}
//				}
//				if (oilCoeff == -1) {
//					constant.put(m, polynomial.coefficient(m));
//				} else {
//					int[] exponents = Arrays.copyOf(m.exponents(), oil + vinegar);
//					exponents[vinegar + oilCoeff] = 0;
//					linear.get(oilCoeff).put(polynomials.getMonomial(exponents), polynomial.coefficient(m));
//				}
//			}
//			List<TExt<C2E>> row = new ArrayList<>();
//			for (Map<Monomial, C2E> coeff : linear) {
//				row.add(ff.getEmbedding(polynomials.getPolynomial(coeff)));
//			}
//			extractedLinear.add(row);
//			extractedConstant.add(ff.getEmbedding(
//					polynomials.subtract(polynomials.getVar(vinegar + i + 1), polynomials.getPolynomial(constant))));
//		}
//		Matrix<TExt<C2E>> matrix = new Matrix<>(extractedLinear);
//		Vector<TExt<C2E>> constant = new Vector<>(extractedConstant);
//		MatrixModule<TExt<C2E>> mm = matrix.getModule(ff);
//		System.out.println(mm.solve(matrix, constant));
//		PolynomialRing<C2E> r = polynomials.eliminateVariable(oil + vinegar);
//		DifferentialForms<C2E> differentialForms = r.differentialForms();
//		PolynomialIdeal<C2E> ideal;
//		do {
//			counter++;
//			if (counter % 10 == 0) {
//				System.out.println("Attempt " + counter);
//			}
////			// byte[] msg = "This is a test message".getBytes(Charset.forName("UTF8"));
//			UOVPrivateKey privateKey = uov.createPrivateKey();
//			Vector<Polynomial<C2E>> publicKey = uov.createPublicKey(privateKey);
//			List<Vector<Polynomial<C2E>>> cotangentRelations = new ArrayList<>();
//			List<Polynomial<C2E>> singularLocus = new ArrayList<>();
//			for (Polynomial<C2E> generator : publicKey.asList()) {
//				Polynomial<C2E> dehomogenized = r.getEmbedding(polynomials.dehomogenize(generator, oil + vinegar));
//				cotangentRelations
//						.add(differentialForms.asGradedVector(differentialForms.derivative(dehomogenized), 1));
//				singularLocus.add(dehomogenized);
//			}
//		Polynomial<C2E> resultant =	r.resultant(singularLocus.get(0), singularLocus.get(1), oil+vinegar-1);
//			System.out.println(resultant);
//			System.out.println(r.eliminateVariable().uniqueFactorization(resultant));
//			Matrix<Polynomial<C2E>> cotangentMatrix = Matrix.fromColumns(cotangentRelations);
//			System.out.println(cotangentMatrix);
//			int rank = (oil + vinegar - 1) - (oil + vinegar - 1 - oil);
//			singularLocus.addAll(cotangentMatrix.getModule(r).minors(cotangentMatrix, rank));
//			ideal = r.getIdeal(singularLocus);
//			GenericProjectiveScheme<C2E> scheme = new GenericProjectiveScheme<>(uov.getField(), polynomials,
//					publicKey.asList());
//			System.out.println(scheme.dimension());
//			System.out.println(scheme.degree());
//			System.out.println(scheme.singularPoints());
		// ideal = r.getZeroIdeal();
//		} while (counter < 1000 && ideal.contains(r.one()));
//		System.out.println(ideal);

//		
//		Polynomial<C2E> p1 = publicKey.get(1);
//		Polynomial<C2E> p2 = publicKey.get(2);
//		Char2Field field = uov.getField();
//		PolynomialRing<C2E> polynomials = uov.getPolynomialRing();
//		int evals = 3;
//		List<C2E> eval = new ArrayList<>();
//		for (int i = 0; i < evals; i++) {
//			eval.add(null);
//		}
//		for (int i = 0; i < polynomials.numberOfVariables() - evals; i++) {
//			eval.add(field.zero()/*getRandomElement()*/);
//		}
//		PolynomialRing<C2E> var3 = AbstractPolynomialRing.getPolynomialRing(field, evals, Monomial.GREVLEX);
//		Polynomial<C2E> v1 = var3.getEmbedding(polynomials.partiallyEvaluate(p1, eval));
//		Polynomial<C2E> v2 = var3.getEmbedding(polynomials.partiallyEvaluate(p2, eval));
//		PolynomialIdeal<C2E> ideal = var3.getIdeal(v1, v2);
//		System.out.println(ideal);
//		System.out.println(ideal.dimension());
//		System.out.println(ideal.degree());
//
//		Polynomial<C2E> r1 = var3.getRandomElement(2);
//		Polynomial<C2E> r2 = var3.getRandomElement(2);
//		PolynomialIdeal<C2E> randomIdeal = var3.getIdeal(r1, r2);
//		System.out.println(randomIdeal);
//		System.out.println(randomIdeal.dimension());
//		System.out.println(randomIdeal.degree());
//
//		CoordinateRing<C2E> cr = ideal.divideOut();
//		AffineScheme<C2E> scheme = new AffineScheme<>(field, cr);
//		System.out.println(scheme);
//		System.out.println(scheme.singularLocus());
//		System.out.println(scheme.singularPoints());
//		System.out.println(scheme.simplify().getSimplification().getRange());
//		AffineNormalization<C2E> normalization = scheme.normalization();
//		System.out.println(normalization.getNormalization());
//		System.out.println(normalization.getNormalization().dimension());
//		System.out.println(normalization.getNormalization().degree());
//		System.out.println(normalization.getNormalization().simplify().getSimplification().getRange());
//
//		CoordinateRing<C2E> randomCr = randomIdeal.divideOut();
//		AffineScheme<C2E> randomScheme = new AffineScheme<>(field, randomCr);
//		System.out.println(randomScheme);
//		System.out.println(randomScheme.singularLocus());
//		System.out.println(randomScheme.singularPoints());
//		System.out.println(randomScheme.simplify().getSimplification().getRange());
//	

//		Polynomial<C2E> s1 = var3.multiply(var3.getVar(1), var3.getVar(2));
//		Polynomial<C2E> s2 = var3.getVar(3);
//		PolynomialIdeal<C2E> compIdeal = var3.getIdeal(s1, s2);
//		System.out.println(compIdeal);
//		System.out.println(compIdeal.dimension());
//		System.out.println(compIdeal.degree());
//		CoordinateRing<C2E> compCr = compIdeal.divideOut();
//		AffineScheme<C2E> compScheme = new AffineScheme<>(field, compCr);
//		System.out.println(compScheme.singularLocus());
//		System.out.println(compScheme.singularPoints());
//		System.out.println(compScheme.isReduced());
//		System.out.println(compScheme.isIrreducible());
//		System.out.println(compScheme.irreducibleComponents());

	}

	@Test
	void randomTest() {
		PrimeField fp = PrimeField.getPrimeField(257);
		PolynomialRing<PFE> var3 = AbstractPolynomialRing.getPolynomialRing(fp, 2, Monomial.GREVLEX);

		Polynomial<PFE> r1 = var3.getRandomElement(2);
		// Polynomial<PFE> r2 = var3.getRandomElement(2);
		PolynomialIdeal<PFE> randomIdeal = var3.getIdeal(r1);
		System.out.println(randomIdeal);
		System.out.println(randomIdeal.dimension());
		System.out.println(randomIdeal.degree());

		CoordinateRing<PFE> randomCr = randomIdeal.divideOut();
		AffineScheme<PFE> randomScheme = new AffineScheme<>(fp, randomCr);
		System.out.println(randomScheme);
		System.out.println(randomScheme.singularLocus());
		System.out.println(randomScheme.singularPoints());
		System.out.println(randomScheme.isReduced());
		System.out.println(randomScheme.isIrreducible());
		System.out.println(randomScheme.irreducibleComponents());
	}

	@Test
	void recoverPublicKeyTest() {
		int oil = 4;
		int vinegar = 10;
		UnbalancedOilAndVinegar uov = new UnbalancedOilAndVinegar(oil, vinegar, 3, Sha3.SHAKE_256);
		UOVPrivateKey privateKey = uov.createPrivateKey();
		Vector<Polynomial<C2E>> publicKey = uov.createPublicKey(privateKey);
		Char2Field field = uov.getField();
		List<Pair<Vector<C2E>, Vector<C2E>>> signatures = new ArrayList<>();
		SecureRandom random = new SecureRandom();
		for (int i = 0; i < (oil + vinegar) * (oil + vinegar); i++) {
			byte[] msg = new byte[32];
			random.nextBytes(msg);
			Vector<C2E> hash = uov.parse(msg);
			Vector<C2E> signature = uov.sign(msg, privateKey);
			assertTrue(uov.verify(msg, signature, publicKey));
			signatures.add(new Pair<>(hash, signature));
		}
		PolynomialRing<C2E> polynomials = AbstractPolynomialRing.getPolynomialRing(field, oil + vinegar,
				Monomial.GREVLEX);
		List<Polynomial<C2E>> basis = new ArrayList<>();
		for (int i = 0; i < oil + vinegar; i++) {
			for (int j = i; j < oil + vinegar; j++) {
				basis.add(polynomials.multiply(polynomials.getVar(i + 1), polynomials.getVar(j + 1)));
			}
		}
		List<Vector<C2E>> images = new ArrayList<>();
		for (Polynomial<C2E> basisPolynomial : basis) {
			List<C2E> image = new ArrayList<>();
			for (Pair<Vector<C2E>, Vector<C2E>> signature : signatures) {
				image.add(polynomials.evaluate(basisPolynomial, signature.getSecond()));
			}
			images.add(new Vector<>(image));
		}
		Matrix<C2E> basisMatrix = Matrix.fromColumns(images);
		MatrixModule<C2E> mm = basisMatrix.getModule(field);
		List<List<C2E>> rightHandSides = new ArrayList<>();
		for (int i = 0; i < oil; i++) {
			rightHandSides.add(new ArrayList<>());
		}
		for (Pair<Vector<C2E>, Vector<C2E>> signature : signatures) {
			for (int i = 0; i < oil; i++) {
				rightHandSides.get(i).add(signature.getFirst().get(i + 1));
			}
		}
		List<Polynomial<C2E>> reconstructed = new ArrayList<>();
		for (int i = 0; i < oil; i++) {
			assertTrue(mm.hasSolution(basisMatrix, new Vector<>(rightHandSides.get(i))));
			Vector<C2E> solution = mm.solve(basisMatrix, new Vector<>(rightHandSides.get(i)));
			Polynomial<C2E> reconstructedPolynomial = polynomials.zero();
			for (int j = 0; j < basis.size(); j++) {
				reconstructedPolynomial = polynomials.add(polynomials.multiply(solution.get(j + 1), basis.get(j)),
						reconstructedPolynomial);
			}
			reconstructed.add(reconstructedPolynomial);
		}
		System.out.println(publicKey);
		System.out.println();
		System.out.println(reconstructed);
		for (Pair<Vector<C2E>, Vector<C2E>> signature : signatures) {
			for (int i = 0; i < oil; i++) {
				assertEquals(signature.getFirst().get(i + 1),
						polynomials.evaluate(reconstructed.get(i), signature.getSecond()));
			}
		}
		for (int i = 0; i < 10; i++) {
			byte[] msg = new byte[32];
			random.nextBytes(msg);
			Vector<C2E> signature = uov.sign(msg, privateKey);
			assertTrue(uov.verify(msg, signature, publicKey));
			assertTrue(uov.verify(msg, signature, new Vector<>(reconstructed)));
		}
	}

}
