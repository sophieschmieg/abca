package fields.numberfields;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import fields.exceptions.InfinityException;
import fields.finitefields.FiniteField.FFE;
import fields.floatingpoint.Complex;
import fields.floatingpoint.Complex.ComplexNumber;
import fields.floatingpoint.Reals;
import fields.helper.AbstractElement;
import fields.helper.AbstractField;
import fields.helper.ExtensionField;
import fields.helper.ExtensionField.ExtensionFieldElement;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.numberfields.NumberField.NFE;
import fields.padics.PAdicExtensionField;
import fields.padics.PAdicExtensionField.PAdicNumberExt;
import fields.padics.PAdicField;
import fields.padics.PAdicField.PAdicNumber;
import fields.polynomials.Monomial;
import fields.polynomials.UnivariatePolynomialRing;
import util.MiscAlgorithms;
import util.Pair;

public class NumberField extends AbstractField<NFE> {
	private Rationals q;
	private ExtensionField<Fraction> field;
	private Reals r;
	private Complex c;
	private List<EmbeddedNumberField> realEmbeddings;
	private List<EmbeddedNumberField> complexEmbeddings;

	public static class NFE extends AbstractElement<NFE> {
		private ExtensionFieldElement<Fraction> e;

		private NFE(ExtensionFieldElement<Fraction> e) {
			this.e = e;
		}

		@Override
		public int compareTo(NFE o) {
			return e.compareTo(o.e);
		}

		public ExtensionFieldElement<Fraction> asExtensionFieldElement() {
			return e;
		}

		public static NFE fromExtensionFieldElement(ExtensionFieldElement<Fraction> t) {
			return new NFE(t);
		}

		@Override
		public String toString() {
			return e.toString();
		}
	}

	public NumberField(Polynomial<Fraction> minimalPolynomial) {
		this.q = Rationals.q();
		this.field = new ExtensionField<>(minimalPolynomial, q);
		this.realEmbeddings = new ArrayList<>();
		this.complexEmbeddings = new ArrayList<>();
		this.r = Reals.r();
		this.c = Complex.c();
		Polynomial<ComplexNumber> minPoly = MiscAlgorithms.mapPolynomial(minimalPolynomial, new MathMap<>() {
			@Override
			public ComplexNumber evaluate(Fraction t) {
				return c.getEmbedding(r.getEmbedding(t));
			}
		}, c.getUnivariatePolynomialRing());
		List<ComplexNumber> roots = c.roots(minPoly);
		for (ComplexNumber alpha : roots) {
			if (alpha.complexPart().equals(r.zero())) {
				this.realEmbeddings.add(new EmbeddedNumberField(this, c.getEmbedding(alpha.realPart())));
			} else if (alpha.complexPart().doubleValue() > 0.0) {
				this.complexEmbeddings.add(new EmbeddedNumberField(this, alpha));
			}
		}
	}

	public NFE getEmbedding(Fraction t) {
		return NFE.fromExtensionFieldElement(field.getEmbedding(t));
	}

	public NFE getEmbedding(IntE t) {
		return getEmbedding(q.getEmbedding(t));
	}

	public NFE getEmbedding(BigInteger t) {
		return getEmbedding(Integers.z().getInteger(t));
	}

	@Override
	public NFE getInteger(BigInteger t) {
		return getEmbedding(t);
	}

	public Polynomial<Fraction> minimalPolynomial() {
		return field.minimalPolynomial();
	}

	public ExtensionField<Fraction> asExtensionField() {
		return field;
	}

	@Override
	public NFE zero() {
		return NFE.fromExtensionFieldElement(field.zero());
	}

	@Override
	public NFE one() {
		return NFE.fromExtensionFieldElement(field.one());
	}

	public NFE alpha() {
		return NFE.fromExtensionFieldElement(field.alpha());
	}

	public boolean isInteger(NFE t) {
		Polynomial<Fraction> minPoly = q.getUnivariatePolynomialRing().normalize(field.minimalPolynomial(t.e));
		for (Monomial m : minPoly.monomials()) {
			if (!q.isInteger(minPoly.coefficient(m))) {
				return false;
			}
		}
		return true;
	}

	@Override
	public BigInteger characteristic() {
		return BigInteger.ZERO;
	}

	@Override
	public NFE add(NFE t1, NFE t2) {
		return NFE.fromExtensionFieldElement(field.add(t1.e, t2.e));
	}

	@Override
	public NFE negative(NFE t) {
		return NFE.fromExtensionFieldElement(field.negative(t.e));
	}

	@Override
	public NFE multiply(NFE t1, NFE t2) {
		return NFE.fromExtensionFieldElement(field.multiply(t1.e, t2.e));
	}

	@Override
	public NFE inverse(NFE t) {
		return NFE.fromExtensionFieldElement(field.inverse(t.e));
	}

	@Override
	public NFE getRandomElement() {
		return NFE.fromExtensionFieldElement(field.getRandomElement());
	}

	@Override
	public boolean isFinite() {
		return false;
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		throw new InfinityException();
	}

	@Override
	public Iterator<NFE> iterator() {
		return new Iterator<NFE>() {
			Iterator<ExtensionFieldElement<Fraction>> it = field.iterator();

			@Override
			public boolean hasNext() {
				return it.hasNext();
			}

			@Override
			public NFE next() {
				return NFE.fromExtensionFieldElement(it.next());
			}
		};
	}

	public List<EmbeddedNumberField> allEmbeddings() {
		List<EmbeddedNumberField> embeddings = new ArrayList<>();
		embeddings.addAll(realEmbeddings());
		embeddings.addAll(complexEmbeddings());
		return embeddings;
	}

	public List<EmbeddedNumberField> realEmbeddings() {
		return Collections.unmodifiableList(this.realEmbeddings);
	}

	public List<EmbeddedNumberField> complexEmbeddings() {
		return Collections.unmodifiableList(this.complexEmbeddings);
	}

	public List<Polynomial<NFE>> factorization(Polynomial<NFE> t) {
		if (t.numberOfVariables() != 1) {
			throw new RuntimeException("Multivariate factorization not implemented.");
		}
		UnivariatePolynomialRing<NFE> ring = this.getUnivariatePolynomialRing();
		List<Polynomial<NFE>> result = new ArrayList<>();
		if (!t.leadingCoefficient().equals(one())) {
			result.add(ring.getEmbedding(t.leadingCoefficient()));
			t = ring.normalize(t);
		}
		List<Polynomial<NFE>> squareFreeFactors = ring.squareFreeFactorization(t);
		for (Polynomial<NFE> sff : squareFreeFactors) {
			result.addAll(factorizeSquareFree(sff));
		}
		return result;
	}

	private List<Polynomial<NFE>> factorizeSquareFree(Polynomial<NFE> t) {
		PolynomialRing<NFE> ring = getUnivariatePolynomialRing();
		Integers z = Integers.z();
		IntE denominator = z.one();
		for (Monomial m : t.monomials()) {
			NFE e = t.coefficient(m);
			for (Fraction coeff : field.asVector(e.e).asList()) {
				IntE c = coeff.getDenominator();
				denominator = z.lcm(denominator, c);
			}
		}
		t = ring.multiply(getEmbedding(denominator), t);
		BigInteger limit = BigInteger.ZERO;
		for (Monomial m : t.monomials()) {
			NFE e = t.coefficient(m);
			for (Fraction coeff : field.asVector(e.e).asList()) {
				IntE c = coeff.getNumerator();
				BigInteger l = c.getValue().abs();
				if (l.compareTo(limit) > 0) {
					limit = l;
				}
			}
		}
		for (Monomial m : this.minimalPolynomial().monomials()) {
			Fraction coeff = this.minimalPolynomial().coefficient(m);
			IntE c = coeff.fraction().getNumerator();
			BigInteger l = c.getValue().abs();
			if (l.compareTo(limit) > 0) {
				limit = l;
			}
		}
		limit = limit.add(BigInteger.ONE);
		limit = limit.shiftLeft(1);
		BigInteger prime = BigInteger.ONE;
		while (true) {
			prime = prime.nextProbablePrime();
			if (BigInteger.valueOf(t.degree()).mod(prime).equals(BigInteger.ZERO)) {
				continue;
			}
			int requiredAccuracy = (int) (Math.log(limit.doubleValue()) / Math.log(prime.doubleValue())) + 1;
			PAdicField base = new PAdicField(prime, 4 * requiredAccuracy);
			Polynomial<PAdicNumber> minPoly = base.fromRational(minimalPolynomial());
			if (!PAdicExtensionField.hasIrreducibleGoodReduction(minPoly, base)) {
				continue;
			}
			PAdicExtensionField zp = new PAdicExtensionField(minPoly, base);
			PolynomialRing<PAdicNumberExt> zpr = zp.getUnivariatePolynomialRing();
			Polynomial<PAdicNumberExt> fp = zp.fromNumberField(this, t);
			Polynomial<FFE> reduced = zp.reduce(fp);
			if (reduced.degree() != fp.degree()) {
				continue;
			}
			PolynomialRing<FFE> reducedRing = zp.reduction().getUnivariatePolynomialRing();
			if (reducedRing.gcd(reduced, reducedRing.derivative(reduced, 1)).degree() > 0) {
				continue;
			}
			System.err.println("Preparation successful");
			List<Polynomial<FFE>> factors = zp.reduction().factorization(reduced);
			System.err.println("Finite Field factorization successful");
			List<Polynomial<PAdicNumberExt>> liftedFactors = new ArrayList<>();
			for (Polynomial<FFE> factor : factors) {
				liftedFactors.add(zp.henselLift(fp, factor, requiredAccuracy));
			}
			System.err.println("Lifts computed successful");

			Map<Integer, Polynomial<PAdicNumberExt>> padicFactors = new TreeMap<>();
			for (int i = 0; i < liftedFactors.size(); i++) {
				padicFactors.put(i, liftedFactors.get(i));
			}
			List<Pair<SortedSet<Integer>, Polynomial<PAdicNumberExt>>> padicCombinedFactors = new ArrayList<>();
			padicCombinedFactors.add(new Pair<>(Collections.emptySortedSet(), zpr.one()));
			List<Polynomial<NFE>> nfeFactors = new ArrayList<>();
			do {
				for (int i = 0; i < 3; i++) {
					if (padicFactors.size() == 0) {
						return nfeFactors;
					}
					CheckCombinationsResult result = checkCombinations(t, zp, requiredAccuracy, padicFactors,
							padicCombinedFactors);
					t = result.cofactor;
					nfeFactors.addAll(result.factors);
					for (int k : result.usedFactors) {
						padicFactors.remove(k);
					}
					padicCombinedFactors = result.combined;
				}
			} while (padicFactors.size() < 3);
			System.err.println("Fast check not successful");

			// for (int i = 0; i< liftedFactors.size();i++) {
			// List<Polynomial<Fraction<IntegerRingElement>>> factors =
			// }

			int gray = 0;
			Polynomial<PAdicNumberExt> testFactor = zpr.one();
			for (int i = 0; i < 1 << liftedFactors.size(); i++) {
				int j = i + 1;
				if (j <= 0) {
					throw new ArithmeticException("too many factors");
				}
				int bit = 0;
				while (j % 2 == 0) {
					bit++;
					j >>= 1;
				}
				if (bit >= liftedFactors.size()) {
					continue;
				}
				if ((gray & (1 << bit)) == 0) {
					testFactor = zpr.multiply(testFactor, liftedFactors.get(bit));
				} else {
					testFactor = zpr.quotientAndRemainder(testFactor, liftedFactors.get(bit)).get(0);
				}
				gray ^= 1 << bit;
				if (testFactor.degree() <= 0 || testFactor.degree() >= t.degree()) {
					continue;
				}
				Polynomial<NFE> factor = MiscAlgorithms.mapPolynomial(testFactor, new MathMap<>() {
					@Override
					public NFE evaluate(PAdicNumberExt number) {
						return zp.roundToNumberField(NumberField.this, number, requiredAccuracy);
					}
				}, ring);
				List<Polynomial<NFE>> qr = ring.quotientAndRemainder(t, factor);
				if (qr.get(1).equals(ring.zero())) {
					List<Polynomial<NFE>> result = new ArrayList<>();
					result.addAll(factorizeSquareFree(factor));
					result.addAll(factorizeSquareFree(qr.get(0)));
					return result;
				}
			}
			return Collections.singletonList(t);
		}
	}

	private static class CheckCombinationsResult {
		private List<Polynomial<NFE>> factors = new ArrayList<>();
		private Polynomial<NFE> cofactor;
		private Set<Integer> usedFactors = new TreeSet<>();
		private List<Pair<SortedSet<Integer>, Polynomial<PAdicNumberExt>>> combined = new ArrayList<>();
	}

	private CheckCombinationsResult checkCombinations(Polynomial<NFE> t, PAdicExtensionField padics,
			int requiredAccuracy, Map<Integer, Polynomial<PAdicNumberExt>> padicFactors,
			List<Pair<SortedSet<Integer>, Polynomial<PAdicNumberExt>>> padicCombinedFactors) {
		CheckCombinationsResult result = new CheckCombinationsResult();
		PolynomialRing<PAdicNumberExt> ring = padics.getUnivariatePolynomialRing();
		result.cofactor = t;
		for (int i : padicFactors.keySet()) {
			Polynomial<PAdicNumberExt> padicFactor = padicFactors.get(i);
			for (Pair<SortedSet<Integer>, Polynomial<PAdicNumberExt>> padicCombinedFactor : padicCombinedFactors) {
				if (padicCombinedFactor.getFirst().size() != 0 && padicCombinedFactor.getFirst().first() <= i) {
					continue;
				}
				SortedSet<Integer> indeces = new TreeSet<>();
				indeces.addAll(padicCombinedFactor.getFirst());
				indeces.add(i);
				Polynomial<PAdicNumberExt> newCombined = ring.multiply(padicFactor, padicCombinedFactor.getSecond());
				CheckFactorResult cfr = checkFactor(result.cofactor, newCombined, padics, requiredAccuracy);
				if (cfr.foundFactor) {
					result.factors.add(cfr.factor);
					result.cofactor = cfr.cofactor;
					result.usedFactors.addAll(indeces);
					break;
				} else {
					result.combined.add(new Pair<>(indeces, newCombined));
				}
			}
		}
		return result;
	}

	private static class CheckFactorResult {
		private boolean foundFactor = false;
		private Polynomial<NFE> factor = null;
		private Polynomial<NFE> cofactor = null;
	}

	private CheckFactorResult checkFactor(Polynomial<NFE> t, Polynomial<PAdicNumberExt> potentialFactor,
			PAdicExtensionField zp, int requiredAccuracy) {
		CheckFactorResult result = new CheckFactorResult();
		PolynomialRing<NFE> ring = getUnivariatePolynomialRing();
		Polynomial<NFE> factor = MiscAlgorithms.mapPolynomial(potentialFactor, new MathMap<>() {
			@Override
			public NFE evaluate(PAdicNumberExt number) {
				return zp.roundToNumberField(NumberField.this, number, requiredAccuracy);
			}
		}, ring);
		List<Polynomial<NFE>> qr = ring.quotientAndRemainder(t, factor);
		if (qr.get(1).equals(ring.zero())) {
			result.foundFactor = true;
			result.factor = factor;
			result.cofactor = qr.get(0);
		}
		return result;
	}
	
	@Override
	public String toString() {
		return field.toString();
	}
}
