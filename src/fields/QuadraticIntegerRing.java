package fields;

import java.math.BigInteger;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import fields.ExtensionField.ExtensionFieldElement;
import fields.FieldOfFractions.Fraction;
import fields.IntegerRing.IntegerRingElement;
import fields.QuadraticIntegerRing.QuadraticInteger;
import util.MiscAlgorithms;

public class QuadraticIntegerRing extends AbstractAlgebra<IntegerRingElement, QuadraticInteger> {
	private BigInteger d;
	private IntegerRing z;
	private FieldOfFractions<IntegerRingElement> q;
	private PolynomialRing<Fraction<IntegerRingElement>> poly;
	private Polynomial<Fraction<IntegerRingElement>> minPoly;
	private ExtensionField<Fraction<IntegerRingElement>> field;
	private QuadraticInteger zero;
	private QuadraticInteger one;
	private QuadraticInteger omegaElement;
	private ExtensionFieldElement<Fraction<IntegerRingElement>> omega;
	private MatrixAlgebra<IntegerRingElement> matrixAlgebra;
	private FreeModule<IntegerRingElement> free;

	public QuadraticIntegerRing(BigInteger d) {
		this.d = d;
		Map<BigInteger, Integer> decomp = MiscAlgorithms.naivePrimeDecomposition(d, d);
		for (int power : decomp.values()) {
			if (power > 1) {
				throw new ArithmeticException("Not square free!");
			}
		}
		z = new IntegerRing();
		matrixAlgebra = new MatrixAlgebra<>(z, 2);
		free = new FreeModule<>(z, 2);
		q = new FieldOfFractions<>(z);
		poly = new PolynomialRing<Fraction<IntegerRingElement>>(q, 1, Polynomial.LEX);
		Polynomial<Fraction<IntegerRingElement>> x2 = poly.getVar(1, 2);
		Polynomial<Fraction<IntegerRingElement>> dp = poly.getEmbedding(q.getEmbedding(z.getInteger(d)));
		minPoly = poly.subtract(x2, dp);
		field = new ExtensionField<>(minPoly, q);
		zero = getEmbedding(z.zero());
		one = getEmbedding(z.one());
		omegaElement = new QuadraticInteger(BigInteger.ZERO, BigInteger.ONE);
		if (d.mod(BigInteger.valueOf(4)).equals(BigInteger.ONE)) {
			omega = field.divide(field.add(field.one(), field.alpha()), field.getInteger(2));
		} else {
			omega = field.alpha();
		}
	}

	@Override
	public Ring<IntegerRingElement> getRing() {
		return z;
	}

	@Override
	public QuadraticInteger zero() {
		return zero;
	}

	@Override
	public QuadraticInteger add(QuadraticInteger s1, QuadraticInteger s2) {
		return new QuadraticInteger(s1.real.add(s2.real), s1.img.add(s2.img));
	}

	@Override
	public QuadraticInteger negative(QuadraticInteger s) {
		return new QuadraticInteger(s.img.negate(), s.img.negate());
	}

	@Override
	public boolean isFree() {
		return true;
	}

	@Override
	public QuadraticInteger getRandomElement() {
		return new QuadraticInteger(z.getRandomElement().getValue(), z.getRandomElement().getValue());
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
	public Iterator<QuadraticInteger> iterator() {
		throw new InfinityException();
	}

	@Override
	public QuadraticInteger one() {
		return one;
	}

	public QuadraticInteger omega() {
		return omegaElement;
	}

	@Override
	public BigInteger characteristic() {
		return BigInteger.ZERO;
	}

	@Override
	public QuadraticInteger multiply(QuadraticInteger t1, QuadraticInteger t2) {
		BigInteger one = BigInteger.ONE;
		BigInteger four = BigInteger.valueOf(4);
		if (d.mod(BigInteger.valueOf(4)).equals(BigInteger.ONE)) {
			return new QuadraticInteger(
					t1.real.multiply(t2.real).add(d.subtract(one).divide(four).multiply(t1.img).multiply(t2.img)),
					t1.real.multiply(t2.img).add(t1.img.multiply(t2.real)).add(t1.img.multiply(t2.img)));
		}
		return new QuadraticInteger(t1.real.multiply(t2.real).add(d.multiply(t1.img).multiply(t2.img)),
				t1.real.multiply(t2.img).add(t1.img.multiply(t2.real)));
	}

	@Override
	public boolean isUnit(QuadraticInteger t) {
		return z.isUnit(norm(t));
	}

	@Override
	public QuadraticInteger inverse(QuadraticInteger t) {
		BigInteger detinv = z.inverse(norm(t)).getValue();
		BigInteger one = BigInteger.ONE;
		BigInteger four = BigInteger.valueOf(4);
		if (d.mod(four).equals(one)) {
			return new QuadraticInteger(detinv.multiply(t.real.add(t.img)), detinv.multiply(t.img.negate()));

		}
		return new QuadraticInteger(detinv.multiply(t.real), detinv.multiply(t.img.negate()));
	}

	@Override
	public boolean isIntegral() {
		return true;
	}

	@Override
	public boolean isZeroDivisor(QuadraticInteger t) {
		return t.equals(zero);
	}

	@Override
	public boolean isEuclidean() {
		return d.compareTo(BigInteger.ZERO) < 0;
	}

	@Override
	public boolean isDivisible(QuadraticInteger dividend, QuadraticInteger divisor) {
		return z.isDivisible(norm(dividend), norm(divisor));
	}

	@Override
	public List<QuadraticInteger> quotientAndRemainder(QuadraticInteger dividend, QuadraticInteger divisor) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public BigInteger euclidMeasure(QuadraticInteger t) {
		return norm(t).getValue();
	}
	
	public IntegerRingElement norm(QuadraticInteger t) {
		return matrixAlgebra.determinant(t.asMatrix());
	}
	
	public IntegerRingElement trace(QuadraticInteger t) {
		return matrixAlgebra.trace(t.asMatrix());
	}

	@Override
	public Iterable<QuadraticInteger> getUnits() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public QuadraticInteger getEmbedding(IntegerRingElement t) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<QuadraticInteger> getGenerators() {
		return Collections.singletonList(omega());
	}

	public class QuadraticInteger implements Element {
		private BigInteger real;
		private BigInteger img;

		public QuadraticInteger(ExtensionFieldElement<Fraction<IntegerRingElement>> element) {
			Vector<Fraction<IntegerRingElement>> asVector = field.asVector(element);
			BigInteger one = BigInteger.ONE;
			BigInteger two = BigInteger.TWO;
			BigInteger four = BigInteger.valueOf(4);
			if (d.mod(four).equals(one)) {
				if (asVector.get(1).getDenominator().equals(z.one())
						&& asVector.get(2).getDenominator().equals(z.one())) {
					img = asVector.get(2).getNumerator().getValue().multiply(two);
					real = asVector.get(1).getNumerator().getValue().subtract(img);
				} else if (asVector.get(1).getDenominator().equals(z.getInteger(2))
						&& asVector.get(2).getDenominator().equals(z.getInteger(2))) {
					img = asVector.get(2).getNumerator().getValue();
					real = asVector.get(1).getNumerator().getValue().subtract(img).divide(two);
				} else {
					throw new ArithmeticException("Not integer!");
				}
			} else {
				if (!asVector.get(1).getDenominator().equals(z.one())
						|| !asVector.get(2).getDenominator().equals(z.one())) {
					throw new ArithmeticException("Not integer!");
				}
				real = asVector.get(1).getNumerator().getValue();
				img = asVector.get(2).getNumerator().getValue();
			}
		}

		public QuadraticInteger(BigInteger real, BigInteger img) {
			this.real = real;
			this.img = img;
		}

		public ExtensionFieldElement<Fraction<IntegerRingElement>> asExtensionFieldElement() {
			ExtensionFieldElement<Fraction<IntegerRingElement>> real = field
					.getEmbedding(q.getEmbedding(z.getInteger(this.real)));
			ExtensionFieldElement<Fraction<IntegerRingElement>> img = field
					.getEmbedding(q.getEmbedding(z.getInteger(this.img)));
			return field.add(real, field.multiply(img, omega));
		}

		public Vector<IntegerRingElement> asVector() {
			return new Vector<>(new IntegerRingElement(real), new IntegerRingElement(img));
		}

		public Matrix<IntegerRingElement> asMatrix() {
			IntegerRingElement[][] m = new IntegerRingElement[2][2];
			BigInteger one = BigInteger.ONE;
			BigInteger four = BigInteger.valueOf(4);
			m[0][0] = new IntegerRingElement(real);
			m[1][0] = new IntegerRingElement(img);
			if (d.mod(four).equals(one)) {
				m[0][1] = new IntegerRingElement(d.subtract(one).divide(four).multiply(img));
				m[1][1] = new IntegerRingElement(real.add(img));
			} else {
				m[0][1] = new IntegerRingElement(d.multiply(img));
				m[1][1] = new IntegerRingElement(real);
			}
			return new Matrix<>(m);
		}

		@Override
		public int compareTo(Element o) {
			QuadraticInteger b = (QuadraticInteger) o;
			int cmp = real.compareTo(b.real);
			if (cmp != 0) {
				return cmp;
			}
			return img.compareTo(b.img);
		}

		@Override
		public boolean equals(Object o) {
			QuadraticInteger b = (QuadraticInteger) o;
			return real.equals(b.real) && img.equals(b.img);
		}

		@Override
		public String toString() {
			return asExtensionFieldElement().toString();
		}
	}
}
