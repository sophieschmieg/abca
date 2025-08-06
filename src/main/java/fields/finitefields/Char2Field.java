package fields.finitefields;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Random;
import java.util.SortedMap;
import java.util.TreeMap;

import fields.exceptions.InfinityException;
import fields.finitefields.Char2Field.C2E;
import fields.finitefields.FiniteField.FFE;
import fields.finitefields.PrimeField.PFE;
import fields.helper.AbstractElement;
import fields.helper.AbstractField;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.UnivariatePolynomial;
import fields.vectors.Vector;
import util.ConCatMap;
import util.FunctionMathMap;
import util.MiscAlgorithms;

public class Char2Field extends AbstractField<C2E> {
	public class C2E extends AbstractElement<C2E> {
		private byte[] value;

		private C2E(byte[] value) {
			this.value = value;
		}

		@Override
		public int compareTo(C2E o) {
			return Arrays.compare(value, o.value);
		}

		public boolean bit(int degree) {
			return (value[degree / 8] & (1 << (degree % 8))) != 0;
		}

		@Override
		public String toString() {
			if (hexPrint) {
				StringBuilder build = new StringBuilder();
				for (int i = value.length - 1; i >= 0; i--) {
					build.append(Integer.toHexString((value[i] & 255) / 16));
					build.append(Integer.toHexString((value[i] & 255) % 16));
				}
				return build.toString();
			}
			boolean zeroOne = true;
			List<String> items = new ArrayList<>();
			for (int i = degree - 1; i >= 0; i--) {
				if (bit(i)) {
					if (i != 0) {
						zeroOne = false;
					}
					items.add("x^" + i);
				}
			}
			if (zeroOne) {
				return bit(0) ? "1" : "0";
			}
			return String.join(" + ", items);
		}

	}

	private int degree;
	private int degreeBytes;
	private byte[] mod;
	private byte[][] powers;
	private C2E zero;
	private C2E one;
	private boolean hexPrint;

	public static Char2Field fromLittleEndian(byte[] mod) {
		return fromLittleEndian(mod, true);
	}

	public static Char2Field fromLittleEndian(byte[] mod, boolean hexPrint) {
		int degree = degree(mod);
		if (degree > 0) {
			return fromLittleEndian(degree, mod, hexPrint);
		}
		throw new ArithmeticException("modulus is zero!");
	}

	public static Char2Field fromLittleEndian(int degree, byte[] mod) {
		return fromLittleEndian(degree, mod, true);
	}

	public static Char2Field fromLittleEndian(int degree, byte[] mod, boolean hexPrint) {
		return new Char2Field(degree, mod, hexPrint);
	}

	public static Char2Field fromMinimalPolynomial(UnivariatePolynomial<PFE> minimalPolynomial) {
		return fromMinimalPolynomial(minimalPolynomial, true);
	}

	public static Char2Field fromMinimalPolynomial(UnivariatePolynomial<PFE> minimalPolynomial, boolean hexPrint) {
		return new Char2Field(minimalPolynomial.degree(), fromPolynomial(minimalPolynomial), hexPrint);
	}

	public static Char2Field fromFiniteField(FiniteField field) {
		return fromFiniteField(field, true);
	}

	public static Char2Field fromFiniteField(FiniteField field, boolean hexPrint) {
		if (!field.characteristic().equals(BigInteger.TWO)) {
			throw new ArithmeticException("Not a characteristic 2 field");
		}
		return fromMinimalPolynomial(field.minimalPolynomial(), hexPrint);
	}

	private static byte[] fromPolynomial(UnivariatePolynomial<PFE> polynomial) {
		PrimeField f2 = PrimeField.getPrimeField(2);
		byte[] result = new byte[MiscAlgorithms.DivRoundUp(polynomial.degree() + 1, 8)];
		for (int i = 0; i <= polynomial.degree(); i++) {
			if (polynomial.univariateCoefficient(i).equals(f2.one())) {
				result[i / 8] |= 1 << (i % 8);
			}
		}
		return result;
	}

	private Char2Field(int degree, byte[] mod, boolean hexPrint) {
		this.hexPrint = hexPrint;
		this.degree = degree;
		this.degreeBytes = MiscAlgorithms.DivRoundUp(degree, 8);
		this.mod = mod;
		this.zero = getElement(new byte[degreeBytes]);
		byte[] oneArray = new byte[degreeBytes];
		oneArray[0] = 1;
		this.one = getElement(oneArray);
		this.powers = new byte[2 * degree][degreeBytes];
		this.powers[0][0] = 1;
		for (int i = 1; i < 2 * degree; i++) {
			this.powers[i] = shift(this.powers[i - 1]);
		}
	}

	@Override
	public String toString() {
		StringBuilder build = new StringBuilder();
		build.append("X^" + degree());
		for (int i = degree - 1; i >= 0; i--) {
			if ((mod[i / 8] & (1 << (i % 8))) != 0) {
				if (i != 0) {
					build.append(" + X^" + i);
				} else {
					build.append(" + 1");
				}
			}
		}
		return "F2[X]/(" + build.toString() + ")";
	}

	public int degree() {
		return degree;
	}

	private byte[] shift(byte[] value) {
		boolean carry = (value[MiscAlgorithms.DivRoundDown(degree - 1, 8)] & 1 << ((degree - 1) % 8)) != 0;
		byte[] result = Arrays.copyOf(value, degreeBytes);
		for (int i = 0; i < degreeBytes; i++) {
			byte prev = i > 0 ? value[i - 1] : 0;
			result[i] = (byte) (((value[i] << 1) & 255) | prev >> 7);
			if (carry) {
				result[i] ^= mod[i];
			}
		}
		return result;
	}

	private byte[] shiftNoMod(byte[] value) {
		int degree = degree(value);
		byte[] result = Arrays.copyOf(value, MiscAlgorithms.DivRoundUp(degree + 2, 8));
		for (int i = 0; i < result.length; i++) {
			byte prev = i > 0 ? value[i - 1] : 0;
			byte current = i < value.length ? value[i] : 0;
			result[i] = (byte) (((current << 1) & 255) | ((prev >> 7) & 1));
		}
		return result;
	}

	private byte[] shiftNoMod(byte[] value, int shift) {
		byte[] result = Arrays.copyOf(value, degreeBytes);
		for (int i = 0; i < shift; i++) {
			result = shiftNoMod(result);
		}
		return result;
	}

	private void xor(byte[] result, byte[] value) {
		if (result.length < value.length) {
			throw new ArithmeticException("Result does not fit!");
		}
		for (int i = 0; i < Math.min(result.length, value.length); i++) {
			result[i] ^= value[i];
		}
	}

	private C2E getElement(byte[] value) {
		if (value.length != degreeBytes) {
			value = Arrays.copyOf(value, degreeBytes);
		}
		if (degree % 8 != 0) {
			value[degree / 8] &= (1 << (degree % 8)) - 1;
		}
		return new C2E(value);
	}

	public C2E fromBytes(byte[] value) {
		return getElement(value);
	}

	public C2E fromFFE(FFE t) {
		PrimeField f2 = PrimeField.getPrimeField(2);
		byte[] result = new byte[degreeBytes];
		for (int i = 0; i <= t.asPolynomial().degree(); i++) {
			if (t.asPolynomial().univariateCoefficient(i).equals(f2.one())) {
				result[i / 8] |= 1 << (i % 8);
			}
		}
		return getElement(result);
	}

	public FFE toFFE(C2E t) {
		PrimeField f2 = PrimeField.getPrimeField(2);
		List<PFE> coeffs = new ArrayList<>();
		for (int i = 0; i < degree; i++) {
			if (t.bit(i)) {
				coeffs.add(f2.one());
			} else {
				coeffs.add(f2.zero());
			}
		}
		return asFiniteField().fromPolynomial(f2.getUnivariatePolynomialRing().getPolynomial(coeffs));
	}

	public FiniteField asFiniteField() {
		PrimeField f2 = PrimeField.getPrimeField(2);
		List<PFE> coeffs = new ArrayList<>();
		for (int i = 0; i < degree + 1; i++) {
			if ((mod[MiscAlgorithms.DivRoundDown(i, 8)] & (1 << (i % 8))) != 0) {
				coeffs.add(f2.one());
			} else {
				coeffs.add(f2.zero());
			}
		}
		return FiniteField.getFiniteField(f2.getUnivariatePolynomialRing().getPolynomial(coeffs), f2);
	}

	@Override
	public C2E zero() {
		return zero;
	}

	@Override
	public C2E one() {
		return one;
	}

	@Override
	public BigInteger characteristic() {
		return BigInteger.TWO;
	}

	@Override
	public C2E add(C2E t1, C2E t2) {
		byte[] result = Arrays.copyOf(t1.value, degreeBytes);
		xor(result, t2.value);
		return getElement(result);
	}

	@Override
	public C2E negative(C2E t) {
		return t;
	}

	public C2E elementFromPolynomial(UnivariatePolynomial<PFE> polynomial) {
		return getElement(fromPolynomial(polynomial));
	}

	@Override
	public C2E multiply(C2E t1, C2E t2) {
		byte[] result = new byte[degreeBytes];
		for (int i = 0; i < degree; i++) {
			for (int j = 0; j < degree; j++) {
				if (t1.bit(i) && t2.bit(j)) {
					xor(result, powers[i + j]);
				}
			}
		}
		return getElement(result);
	}

	@Override
	public FactorizationResult<Polynomial<C2E>, C2E> factorization(UnivariatePolynomial<C2E> t) {
		FiniteField fq = asFiniteField();
		FactorizationResult<Polynomial<FFE>, FFE> factorization = fq.factorization(
				fq.getUnivariatePolynomialRing().getEmbedding(t, new FunctionMathMap<>((C2E s) -> toFFE(s))));
		SortedMap<Polynomial<C2E>, Integer> result = new TreeMap<>();
		for (Polynomial<FFE> factor : factorization.primeFactors()) {
			result.put(getUnivariatePolynomialRing().getEmbedding(factor, new FunctionMathMap<>((FFE s) -> fromFFE(s))),
					factorization.multiplicity(factor));
		}
		return new FactorizationResult<>(fromFFE(factorization.getUnit()), result);
	}

	@Override
	public Extension<C2E, PFE, FFE, FiniteField> getExtension(UnivariatePolynomial<C2E> minimalPolynomial) {
		MathMap<C2E, FFE> embedding = new FunctionMathMap<>((C2E t) -> toFFE(t));
		FiniteField asField = asFiniteField();
		Extension<FFE, PFE, FFE, FiniteField> extension = asField
				.getExtension(asField.getUnivariatePolynomialRing().getEmbedding(minimalPolynomial, embedding));
		return new Extension<>(extension.extension(), this, new ConCatMap<>(embedding, extension.embeddingMap()),
				new ConCatMap<>(extension.asVectorMap(), new FunctionMathMap<>(
						(Vector<FFE> t) -> Vector.mapVector(new FunctionMathMap<>((FFE s) -> fromFFE(s)), t))));
	}

	private static int degree(byte[] polynomial) {
		for (int index = polynomial.length - 1; index >= 0; index--) {
			if (polynomial[index] == 0) {
				continue;
			}
			int modDegree = 0;
			int modCopy = polynomial[index] & 255;
			while (modCopy != 1) {
				modCopy >>= 1;
				modDegree++;
			}
			return 8 * index + modDegree;
		}
		return -1;
	}

	private byte[][] quotientAndRemainder(byte[] dividend, byte[] divisor) {
		byte[] mod = dividend;
		int divisorDegree = degree(divisor);
		int modDegree = degree(mod);
		byte[] quotient = new byte[MiscAlgorithms.DivRoundUp(modDegree - divisorDegree + 1, 8)];
		while (modDegree >= divisorDegree) {
			byte[] shifted = shiftNoMod(divisor, modDegree - divisorDegree);
			quotient[MiscAlgorithms.DivRoundDown(modDegree - divisorDegree,
					8)] |= (byte) (1 << ((modDegree - divisorDegree) % 8));
			xor(shifted, mod);
			modDegree = degree(shifted);
			mod = Arrays.copyOf(shifted, MiscAlgorithms.DivRoundUp(modDegree + 1, 8));
		}
		byte[][] result = new byte[2][];
		result[0] = quotient;
		result[1] = mod;
		return result;
	}

	@Override
	public C2E inverse(C2E t) {
		// xi * a + yi * b = ri
		// x0 = 1, y0 = 0
		// x1 = 0, y0 = 1
		// ri-1 = qi * ri + ri+1
		// xi-1*a + yi-1*b = qi * (xi*a + yi*b) + xi+1*a + yi+1*b
		// xi+1 = (xi-1 * qi*xi)
		byte[] a = mod;
		byte[] b = t.value;
		C2E prev = zero();
		C2E current = one();
		while (degree(b) > 0) {
			byte[][] qr = quotientAndRemainder(a, b);
			a = b;
			b = qr[1];
			C2E next = add(prev, multiply(current, getElement(qr[0])));
			prev = current;
			current = next;
		}
		return current;
	}

	@Override
	public Exactness exactness() {
		return Exactness.EXACT;
	}

	@Override
	public C2E getRandomElement() {
		byte[] result = new byte[degreeBytes];
		Random rng = new Random();
		rng.nextBytes(result);
		if (degree % 8 != 0) {
			result[degreeBytes - 1] &= (1 << (degree % 8)) - 1;
		}
		return getElement(result);
	}

	@Override
	public boolean isFinite() {
		return true;
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		return BigInteger.ONE.shiftLeft(degree);
	}

	@Override
	public Iterator<C2E> iterator() {
		return new Iterator<>() {
			private BigInteger next = BigInteger.ZERO;

			@Override
			public C2E next() {
				byte[] result = new byte[degreeBytes];
				for (int i = 0; i < degree; i++) {
					if (next.testBit(i)) {
						result[MiscAlgorithms.DivRoundDown(i, 8)] |= (byte) (1 << (i % 8));
					}
				}
				next = next.add(BigInteger.ONE);
				return getElement(result);
			}

			@Override
			public boolean hasNext() {
				return next.compareTo(getNumberOfElements()) <= 0;
			}
		};
	}

	@Override
	public boolean hasCharacteristicRoot(C2E t, int power) {
		return true;
	}

	@Override
	public C2E characteristicRoot(C2E t, int power) {
		return power(t, characteristic().pow(degree() - (power % degree())));
	}

}
