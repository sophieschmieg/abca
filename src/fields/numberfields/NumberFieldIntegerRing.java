package fields.numberfields;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import fields.exceptions.InfinityException;
import fields.helper.AbstractAlgebra;
import fields.helper.AbstractIdeal;
import fields.helper.ExtensionField;
import fields.helper.ExtensionField.ExtensionFieldElement;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.Element;
import fields.interfaces.Ideal;
import fields.interfaces.Polynomial;
import fields.interfaces.Ring;
import fields.numberfields.NumberFieldIntegerRing.NumberFieldInteger;
import fields.polynomials.Monomial;
import fields.vectors.FreeModule;
import fields.vectors.MatrixAlgebra;
import fields.vectors.Vector;

public class NumberFieldIntegerRing extends AbstractAlgebra<IntE, NumberFieldInteger> {
	private Integers z;
	private Rationals q;
	private ExtensionField<Fraction> field;
	private MatrixAlgebra<IntE> matrixAlgebra;
	private FreeModule<IntE> free;

	public NumberFieldIntegerRing(ExtensionField<Fraction> field) {
		this.z = Integers.z();
		this.q = Rationals.q();
		this.field = field;
		//this.matrixAlgebra = new MatrixAlgebra<>(this.z, this.field.degree());
		this.free = new FreeModule<>(this.z, this.field.degree());
	}

	@Override
	public Ring<IntE> getRing() {
		return z;
	}

	public int degree() {
		return this.field.degree();
	}

	@Override
	public NumberFieldInteger zero() {
		return getInteger(0);
	}

	@Override
	public NumberFieldInteger add(NumberFieldInteger s1, NumberFieldInteger s2) {
		return getElementNoCheck(field.add(s1.fieldElement, s2.fieldElement));
	}

	@Override
	public NumberFieldInteger negative(NumberFieldInteger s) {
		return getElementNoCheck(field.negative(s.fieldElement));
	}

	@Override
	public boolean isFree() {
		return true;
	}

	@Override
	public NumberFieldInteger getRandomElement() {
		throw new InfinityException();
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
	public Iterator<NumberFieldInteger> iterator() {
		throw new InfinityException();
	}

	@Override
	public NumberFieldInteger one() {
		return getInteger(1);
	}

	@Override
	public BigInteger characteristic() {
		return BigInteger.ZERO;
	}

	@Override
	public NumberFieldInteger multiply(NumberFieldInteger t1, NumberFieldInteger t2) {
		return getElementNoCheck(field.multiply(t1.fieldElement, t2.fieldElement));
	}

	@Override
	public boolean isUnit(NumberFieldInteger t) {
		return norm(t).equals(BigInteger.ONE);
	}

	@Override
	public NumberFieldInteger inverse(NumberFieldInteger t) {
		return getElement(field.inverse(t.fieldElement));
	}

	@Override
	public boolean isIntegral() {
		return true;
	}

	@Override
	public boolean isZeroDivisor(NumberFieldInteger t) {
		return t.equals(zero());
	}

	@Override
	public boolean isEuclidean() {
		return false;
	}

	public BigInteger norm(NumberFieldInteger t) {
		return this.field.norm(t.fieldElement).getNumerator().getValue();
	}

	public boolean isInteger(ExtensionFieldElement<Fraction> t) {
		Polynomial<Fraction> minPoly = field.minimalPolynomial(t);
		if (!minPoly.leadingCoefficient().equals(q.one())) {
			return false;
		}
		for (int i = 0; i < degree(); i++) {
			Monomial m = q.getUnivariatePolynomialRing().getMonomial(new int[] { i });
			if (!minPoly.coefficient(m).getDenominator().equals(z.one())) {
				return false;
			}
		}
		return true;
	}

	public NumberFieldInteger getElement(ExtensionFieldElement<Fraction> t) {
		if (!isInteger(t)) {
			throw new ArithmeticException("not integer");
		}
		return new NumberFieldInteger(t);
	}

	private NumberFieldInteger getElementNoCheck(ExtensionFieldElement<Fraction> t) {
		return new NumberFieldInteger(t);
	}

	@Override
	public boolean isDivisible(NumberFieldInteger dividend, NumberFieldInteger divisor) {
		return isInteger(field.divide(dividend.fieldElement, divisor.fieldElement));
	}

	@Override
	public List<NumberFieldInteger> quotientAndRemainder(NumberFieldInteger dividend, NumberFieldInteger divisor) {
		if (!isDivisible(dividend, divisor)) {
			throw new ArithmeticException("Not divisible!");
		}
		List<NumberFieldInteger> result = new ArrayList<>();
		result.add(getElementNoCheck(field.divide(dividend.fieldElement, divisor.fieldElement)));
		result.add(zero());
		return result;
	}

	@Override
	public BigInteger euclidMeasure(NumberFieldInteger t) {
		throw new ArithmeticException("Not Euclidean!");
	}
	
	public NumberFieldInteger projectToUnit(NumberFieldInteger t) {
		return one();
	}
	
	@Override
	public Iterable<NumberFieldInteger> getUnits() {
		throw new InfinityException();
	}

	public NumberFieldInteger getInteger(int t) {
		return getEmbedding(z.getInteger(t));
	}

	public NumberFieldInteger getInteger(BigInteger t) {
		return getEmbedding(z.getInteger(t));
	}

	@Override
	public NumberFieldInteger getEmbedding(IntE t) {
		return getElementNoCheck(field.getEmbedding(q.getEmbedding(t)));
	}

	@Override
	public Ideal<NumberFieldInteger> getIdeal(List<NumberFieldInteger> generators) {
		return new NumberFieldIdeal(generators);
	}
	
	public Ideal<NumberFieldInteger> intersect(Ideal<NumberFieldInteger> t1, Ideal<NumberFieldInteger> t2) {
		return multiply(t1, t2);
	}
	
	public int krullDimension() {
		throw new UnsupportedOperationException("Not implemented");
	}

	public static class NumberFieldInteger implements Element<NumberFieldInteger> {
		private ExtensionFieldElement<Fraction> fieldElement;
		private Polynomial<Fraction> minPoly;

		private NumberFieldInteger(ExtensionFieldElement<Fraction> fieldElement) {
			this.fieldElement = fieldElement;
		}

		@Override
		public boolean equals(Object o) {
			NumberFieldInteger other = (NumberFieldInteger) o;
			return this.fieldElement.equals(other.fieldElement);
		}

		@Override
		public int compareTo(NumberFieldInteger o) {
			return this.fieldElement.compareTo(o.fieldElement);
		}
	}

	public class NumberFieldIdeal extends AbstractIdeal<NumberFieldInteger> {
		private List<NumberFieldInteger> generators;
		public NumberFieldIdeal(List<NumberFieldInteger> generators) {
			super (NumberFieldIntegerRing.this);
			this.generators = generators;
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
		public List<NumberFieldInteger> generators() {
			return generators;
		}
		
		@Override
		public List<NumberFieldInteger> generate(NumberFieldInteger t) {
			throw new UnsupportedOperationException("Not implemented");
		}

		@Override
		public NumberFieldInteger residue(NumberFieldInteger t) {
			throw new UnsupportedOperationException("Not implemented");
		}

		@Override
		public boolean isPrime() {
			// TODO Auto-generated method stub
			return false;
		}

		@Override
		public boolean isMaximal() {
			// TODO Auto-generated method stub
			return false;
		}

		@Override
		public boolean contains(NumberFieldInteger t) {
			// TODO Auto-generated method stub
			return false;
		}

	}

	@Override
	public List<NumberFieldInteger> getAlgebraGenerators() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public boolean isLinearIndependent(List<NumberFieldInteger> s) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public boolean isGeneratingModule(List<NumberFieldInteger> s) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public List<NumberFieldInteger> getModuleGenerators() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Vector<IntE> asVector(NumberFieldInteger s) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public boolean isGeneratingAlgebra(List<NumberFieldInteger> s) {
		// TODO Auto-generated method stub
		return false;
	}
}
