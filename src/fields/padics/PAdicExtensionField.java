package fields.padics;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import fields.exceptions.InfinityException;
import fields.finitefields.FiniteField;
import fields.finitefields.FiniteField.FFE;
import fields.finitefields.PrimeField.PrimeFieldElement;
import fields.helper.AbstractElement;
import fields.helper.AbstractField;
import fields.helper.ExtensionField;
import fields.helper.ExtensionField.ExtensionFieldElement;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.ValueField;
import fields.numberfields.NumberField;
import fields.numberfields.NumberField.NFE;
import fields.padics.PAdicExtensionField.PAdicNumberExt;
import fields.padics.PAdicField.PAdicNumber;
import fields.polynomials.AbstractPolynomialRing;
import util.MiscAlgorithms;

public class PAdicExtensionField extends AbstractField<PAdicNumberExt> implements ValueField<PAdicNumberExt> {
	private PAdicField base;
	private ExtensionField<PAdicNumber> field;
	private FiniteField reduction;
	private Polynomial<PAdicNumber> minPoly;

	public static class PAdicNumberExt extends AbstractElement<PAdicNumberExt> {
		private ExtensionFieldElement<PAdicNumber> e;

		private PAdicNumberExt(ExtensionFieldElement<PAdicNumber> e) {
			this.e = e;
		}

		@Override
		public int compareTo(PAdicNumberExt o) {
			return e.compareTo(o.e);
		}

		public ExtensionFieldElement<PAdicNumber> asExtensionFieldElement() {
			return e;
		}

		public static PAdicNumberExt fromExtensionFieldElement(ExtensionFieldElement<PAdicNumber> t) {
			return new PAdicNumberExt(t);
		}

		@Override
		public String toString() {
			return e.toString();
		}

	}

	public PAdicExtensionField(PAdicField base) {
		this(base, 1);
	}

	public PAdicExtensionField(PAdicField base, int degree) {
		this.base = base;
		this.reduction = new FiniteField(this.base.reduction(), degree);
		this.minPoly = this.base.lift(this.reduction.minimalPolynomial());
		this.field = new ExtensionField<PAdicNumber>(this.minPoly, this.base);
	}

	public static boolean hasIrreducibleGoodReduction(Polynomial<PAdicNumber> minimalPolynomial, PAdicField base) {
		Polynomial<PrimeFieldElement> reducedMinimalPolynomial = base.reduce(minimalPolynomial);
		return minimalPolynomial.degree() == reducedMinimalPolynomial.degree()
				&& base.reduction().factorization(reducedMinimalPolynomial).size() == 1;
	}

	public PAdicExtensionField(Polynomial<PAdicNumber> minimalPolynomial, PAdicField base) {
		if (!hasIrreducibleGoodReduction(minimalPolynomial, base)) {
			throw new ArithmeticException("Minimal polynomial does not have good reduction");
		}
		this.base = base;
		this.reduction = new FiniteField(base.reduce(minimalPolynomial), this.base.reduction());
		this.minPoly = minimalPolynomial;
		this.field = new ExtensionField<PAdicNumber>(this.minPoly, this.base);
	}

	public PAdicExtensionField(BigInteger prime, int degree, int maxAccuracy) {
		this(new PAdicField(prime, maxAccuracy), degree);
	}

	public PAdicNumberExt getEmbedding(PAdicNumber t) {
		return new PAdicNumberExt(field.getEmbedding(t));
	}

	public Polynomial<PAdicNumber> minimalPolynomial() {
		return this.minPoly;
	}

	public ExtensionField<PAdicNumber> asExtensionField() {
		return field;
	}

	@Override
	public PAdicNumberExt zero() {
		return new PAdicNumberExt(field.zero());
	}

	@Override
	public PAdicNumberExt one() {
		return new PAdicNumberExt(field.one());
	}

	@Override
	public double value(PAdicNumberExt t) {
		return Math.pow(base.value(field.norm(t.e)), 1.0 / field.degree());
	}

	@Override
	public BigInteger characteristic() {
		return BigInteger.ZERO;
	}

	@Override
	public PAdicNumberExt add(PAdicNumberExt t1, PAdicNumberExt t2) {
		return new PAdicNumberExt(field.add(t1.e, t2.e));
	}

	@Override
	public PAdicNumberExt negative(PAdicNumberExt t) {
		return new PAdicNumberExt(field.negative(t.e));
	}

	@Override
	public PAdicNumberExt multiply(PAdicNumberExt t1, PAdicNumberExt t2) {
		return new PAdicNumberExt(field.multiply(t1.e, t2.e));
	}

	@Override
	public PAdicNumberExt inverse(PAdicNumberExt t) {
		return new PAdicNumberExt(field.inverse(t.e));
	}

	@Override
	public PAdicNumberExt getRandomElement() {
		return new PAdicNumberExt(field.getRandomElement());
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
	public Iterator<PAdicNumberExt> iterator() {
		return new Iterator<PAdicNumberExt>() {
			private Iterator<ExtensionFieldElement<PAdicNumber>> it = field.iterator();

			@Override
			public boolean hasNext() {
				return it.hasNext();
			}

			@Override
			public PAdicNumberExt next() {
				return new PAdicNumberExt(it.next());
			}
		};
	}

	public PAdicNumberExt fromRational(Fraction t) {
		return getEmbedding(base.fromRational(t));
	}

	public Polynomial<PAdicNumberExt> fromRational(Polynomial<Fraction> t) {
		PolynomialRing<Fraction> rationalRing = t.getPolynomialRing();
		PolynomialRing<PAdicNumberExt> ring = AbstractPolynomialRing.getPolynomialRing(this,
				rationalRing.numberOfVariables(), rationalRing.getComparator());
		return MiscAlgorithms.mapPolynomial(t, new MathMap<>() {
			@Override
			public PAdicNumberExt evaluate(Fraction n) {
				return fromRational(n);
			}
		}, ring);
	}

	public boolean numberFieldConsistent(NumberField f) {
		return this.minPoly.equals(base.fromRational(f.minimalPolynomial()));
	}

	public PAdicNumberExt fromNumberField(NumberField f, NFE t) {
		return PAdicNumberExt.fromExtensionFieldElement(field
				.fromPolynomial(base.fromRational(f.asExtensionField().asPolynomial(t.asExtensionFieldElement()))));
	}

	public Polynomial<PAdicNumberExt> fromNumberField(NumberField f, Polynomial<NFE> t) {
		PolynomialRing<NFE> numberFieldRing = t.getPolynomialRing();
		PolynomialRing<PAdicNumberExt> ring = AbstractPolynomialRing.getPolynomialRing(this,
				numberFieldRing.numberOfVariables(), numberFieldRing.getComparator());
		return MiscAlgorithms.mapPolynomial(t, new MathMap<>() {
			@Override
			public PAdicNumberExt evaluate(NFE n) {
				return fromNumberField(f, n);
			}
		}, ring);
	}

	public FiniteField reduction() {
		return reduction;
	}

	public FFE reduce(PAdicNumberExt t) {
		return FFE.fromExtensionFieldElement(
				reduction.asExtensionField().fromPolynomial(base.reduce(field.asPolynomial(t.e))));
	}

	public PAdicNumberExt lift(FFE t) {
		return PAdicNumberExt.fromExtensionFieldElement(field
				.fromPolynomial(base.lift(reduction.asExtensionField().asPolynomial(t.asExtensionFieldElement()))));
	}

	public Polynomial<FFE> reduce(Polynomial<PAdicNumberExt> t) {
		PolynomialRing<PAdicNumberExt> ring = t.getPolynomialRing();
		PolynomialRing<FFE> reducedRing = AbstractPolynomialRing.getPolynomialRing(reduction, ring.numberOfVariables(),
				ring.getComparator());
		return MiscAlgorithms.mapPolynomial(t, new MathMap<>() {
			@Override
			public FFE evaluate(PAdicNumberExt n) {
				return reduce(n);
			}
		}, reducedRing);
	}

	public Polynomial<PAdicNumberExt> lift(Polynomial<FFE> t) {
		PolynomialRing<FFE> ring = t.getPolynomialRing();
		PolynomialRing<PAdicNumberExt> padicRing = AbstractPolynomialRing.getPolynomialRing(this,
				ring.numberOfVariables(), ring.getComparator());
		return MiscAlgorithms.mapPolynomial(t, new MathMap<>() {
			@Override
			public PAdicNumberExt evaluate(FFE n) {
				return lift(n);
			}
		}, padicRing);
	}

	public PAdicNumberExt round(PAdicNumberExt t, int accuracy) {
		return PAdicNumberExt.fromExtensionFieldElement(
				field.fromPolynomial(base.round(field.asPolynomial(t.asExtensionFieldElement()), accuracy)));
	}

	public Polynomial<PAdicNumberExt> round(Polynomial<PAdicNumberExt> t, int accuracy) {
		return MiscAlgorithms.mapPolynomial(t, new MathMap<>() {
			@Override
			public PAdicNumberExt evaluate(PAdicNumberExt n) {
				return round(n, accuracy);
			}
		}, t.getPolynomialRing());
	}

	public BigInteger roundToInteger(PAdicNumberExt t, int accuracy) {
		return base
				.roundToInteger(field.asPolynomial(t.e)
						.coefficient(base.getUnivariatePolynomialRing().getMonomial(new int[] { 0 })), accuracy)
				.getValue();
	}

	public NFE roundToNumberField(NumberField f, PAdicNumberExt t, int accuracy) {
		Polynomial<PAdicNumber> asPoly = field.asPolynomial(t.asExtensionFieldElement());
		Polynomial<Fraction> roundedPoly = MiscAlgorithms.mapPolynomial(asPoly, new MathMap<>() {

			@Override
			public Fraction evaluate(PAdicNumber t) {
				return Rationals.q().getEmbedding(base.roundToInteger(t, accuracy));
			}
		}, Rationals.q().getUnivariatePolynomialRing());
		return NFE.fromExtensionFieldElement(f.asExtensionField().fromPolynomial(roundedPoly));
	}

	public PAdicNumberExt henselLift(Polynomial<PAdicNumberExt> f, FFE aReduced, int accuracy) {
		Polynomial<FFE> fReduced = reduce(f);
		PolynomialRing<FFE> reducedRing = reduction.getUnivariatePolynomialRing();
		PolynomialRing<PAdicNumberExt> ring = getUnivariatePolynomialRing();
		if (fReduced.degree() != f.degree()) {
			throw new ArithmeticException("Leading coefficient divisible by p!");
		}
		if (!reducedRing.evaluate(fReduced, aReduced).equals(reduction.zero())) {
			throw new ArithmeticException("reduced polynomial does not have a zero at a!");
		}
		if (reducedRing.evaluate(reducedRing.derivative(fReduced, 1), aReduced).equals(reduction.zero())) {
			throw new ArithmeticException("reduced derivative does have a zero at a!");
		}
		PAdicNumberExt a = lift(aReduced);
		Polynomial<PAdicNumberExt> derivative = ring.derivative(f, 1);
		int achievedAccuracy = 1;
		while (achievedAccuracy < accuracy) {
			a = subtract(a, divide(ring.evaluate(f, a), ring.evaluate(derivative, a)));
			achievedAccuracy *= 2;
		}
		return a;
	}

	public Polynomial<PAdicNumberExt> henselLift(Polynomial<PAdicNumberExt> f, Polynomial<FFE> gReduced, int accuracy) {
		PolynomialRing<PAdicNumberExt> ring = getUnivariatePolynomialRing();
		Polynomial<FFE> fReduced = reduce(f);
		if (fReduced.degree() != f.degree()) {
			throw new ArithmeticException("Leading coefficient divisible by p!");
		}
		PolynomialRing<FFE> reducedRing = reduction.getUnivariatePolynomialRing();
		List<Polynomial<FFE>> qr = reducedRing.quotientAndRemainder(fReduced, gReduced);
		Polynomial<FFE> hReduced = qr.get(0);
		if (!qr.get(1).equals(reducedRing.zero())) {
			throw new ArithmeticException("Not a divisor!");
		}
		List<Polynomial<FFE>> ee = reducedRing.extendedEuclidean(gReduced, hReduced);
		if (ee.get(0).degree() != 0) {
			throw new ArithmeticException("Not coprime!");
		}
		FFE gcd = reduction.inverse(ee.get(0).leadingCoefficient());
		Polynomial<FFE> aReduced = reducedRing.multiply(gcd, ee.get(1));
		Polynomial<FFE> bReduced = reducedRing.multiply(gcd, ee.get(2));
		int degree = f.degree();
		Polynomial<PAdicNumberExt> g = lift(gReduced);
		Polynomial<PAdicNumberExt> h = ring.divide(f, g);//lift(hReduced);
		List<Polynomial<PAdicNumberExt>> eeLifted = ring.extendedEuclidean(g, h);
		if (ee.get(0).degree() != 0) {
			throw new ArithmeticException("Not coprime!");
		}
		PAdicNumberExt gcdLifted = inverse(eeLifted.get(0).leadingCoefficient());
		Polynomial<PAdicNumberExt> a = ring.multiply(gcdLifted, eeLifted.get(1));//lift(aReduced);
		Polynomial<PAdicNumberExt> b = ring.multiply(gcdLifted, eeLifted.get(2));//lift(bReduced);
		int achievedAccuracy = 1;
		while (achievedAccuracy < 4 * accuracy) {
			achievedAccuracy *= 2;
			Polynomial<PAdicNumberExt> q = ring.round(round(ring.subtract(f, ring.multiply(g, h)), achievedAccuracy), degree);
			Polynomial<PAdicNumberExt> gOld = g;
			g = ring.round(round(ring.add(g, ring.multiply(b, q)), achievedAccuracy), degree);
			//if (achievedAccuracy < 2 * accuracy) {
			Polynomial<PAdicNumberExt> hOld = h;
			h = ring.round(round(ring.add(h, ring.multiply(a, q)), achievedAccuracy), degree);
			Polynomial<PAdicNumberExt> r = ring
						.round(round(ring.subtract(ring.add(ring.multiply(a, g), ring.multiply(b, h)), ring.one()), achievedAccuracy), degree);
				a = ring.round(round(ring.multiply(a, ring.subtract(ring.one(), r)), achievedAccuracy), degree);
				b = ring.round(round(ring.multiply(b, ring.subtract(ring.one(), r)), achievedAccuracy), degree);
			//}
		}
		g = round(g, accuracy);
		return g;
	}

	@Override
	public List<Polynomial<PAdicNumberExt>> factorization(Polynomial<PAdicNumberExt> t) {
		if (!t.leadingCoefficient().equals(one())) {
			throw new ArithmeticException("not normalized");
		}
		if (t.degree() <= 1) {
			return Collections.singletonList(t);
		}
		if (this.base.getMaxAccuracy() == PAdicField.INFINITE_ACCURACY) {
			throw new UnsupportedOperationException("Cannot do arbitrary accuracy factorization over p adics");
		}
		List<Polynomial<FFE>> reducedFactors = reduction().factorization(reduce(t));
		List<Polynomial<PAdicNumberExt>> factors = new ArrayList<>();
		Polynomial<PAdicNumberExt> f = t;
		for (Polynomial<FFE> reducedFactor : reducedFactors) {
			Polynomial<PAdicNumberExt> factor = henselLift(f, reducedFactor, base.getMaxAccuracy());
			factors.add(factor);
			f = getUnivariatePolynomialRing().quotientAndRemainder(f, factor).get(0);
		}
		return factors;
	}

	@Override
	public String toString() {
		return field.toString();
	}
}
