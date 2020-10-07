package fields.numberfields;

import java.math.BigInteger;
import java.util.Iterator;
import java.util.List;

import fields.exceptions.InfinityException;
import fields.floatingpoint.Complex;
import fields.floatingpoint.Complex.ComplexNumber;
import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.helper.AbstractField;
import fields.helper.ExtensionField;
import fields.integers.Integers.IntE;
import fields.integers.Rationals.Fraction;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.ValueField;
import fields.numberfields.NumberField.NFE;
import util.MiscAlgorithms;

public class EmbeddedNumberField extends AbstractField<NFE> implements ValueField<NFE> {
	private NumberField field;
	private Reals r;
	private Complex c;
	private ComplexNumber alpha;
	private boolean realEmbedding;
	
	EmbeddedNumberField(NumberField field, ComplexNumber alpha) {
		this.field = field;
		this.r = Reals.r();
		this.c = Complex.c();
		this.alpha = alpha;
		this.realEmbedding = alpha.complexPart().equals(r.zero());
	}

	public NFE getEmbedding(Fraction t) {
		return field.getEmbedding(t);
	}

	public NFE getEmbedding(IntE t) {
		return field.getEmbedding(t);
	}

	public NFE getEmbedding(BigInteger t) {
		return field.getEmbedding(t);
	}

	@Override
	public NFE getInteger(BigInteger t) {
		return field.getInteger(t);
	}

	public Polynomial<Fraction> minimalPolynomial() {
		return field.minimalPolynomial();
	}

	public ExtensionField<Fraction> asExtensionField() {
		return field.asExtensionField();
	}

	@Override
	public NFE zero() {
		return field.zero();
	}

	@Override
	public NFE one() {
		return field.one();
	}

	public NFE alpha() {
		return field.alpha();
	}

	public boolean isInteger(NFE t) {
		return field.isInteger(t);
	}

	@Override
	public NFE add(NFE t1, NFE t2) {
		return field.add(t1, t2);
	}

	@Override
	public NFE negative(NFE t) {
		return field.negative(t);
	}

	@Override
	public NFE multiply(NFE t1, NFE t2) {
		return field.multiply(t1, t2);
	}

	@Override
	public NFE inverse(NFE t) {
		return field.inverse(t);
	}

	@Override
	public NFE getRandomElement() {
		return field.getRandomElement();
	}

	@Override
	public Iterator<NFE> iterator() {
		return field.iterator();
	}

	@Override
	public List<Polynomial<NFE>> factorization(Polynomial<NFE> t) {
		return field.factorization(t);
	}

	@Override
	public BigInteger characteristic() {
		return BigInteger.ZERO;
	}

	@Override
	public boolean isFinite() {
		return false;
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		throw new InfinityException();
	}

	public boolean isRealEmbedding() {
		return realEmbedding;
	}
	
	public Real asRealNumber(NFE t) {
		if (!isRealEmbedding()) {
			throw new ArithmeticException("Not a real embedding");
		}
		return asComplexNumber(t).realPart();
	}
	
	public ComplexNumber asComplexNumber(NFE t) {
		Polynomial<ComplexNumber> asPolynomial = MiscAlgorithms
				.mapPolynomial(field.asExtensionField().asPolynomial(t.asExtensionFieldElement()), new MathMap<>() {
					@Override
					public ComplexNumber evaluate(Fraction q) {
						return c.getEmbedding(q);
					}
				}, c.getUnivariatePolynomialRing());
		return c.getUnivariatePolynomialRing().evaluate(asPolynomial, alpha);
	}
	
	@Override
	public double value(NFE t) {
		return c.value(asComplexNumber(t));
	}
	

	@Override
	public String toString() {
		return field.toString() + ", alpha=" + alpha;
	}
}
