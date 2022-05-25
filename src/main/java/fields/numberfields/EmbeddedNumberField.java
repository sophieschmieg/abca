package fields.numberfields;

import java.math.BigInteger;
import java.util.Collections;
import java.util.Iterator;

import fields.exceptions.InfinityException;
import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.helper.AbstractField;
import fields.integers.Integers.IntE;
import fields.integers.Rationals.Fraction;
import fields.interfaces.Element;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.ValueField;
import fields.numberfields.NumberField.NFE;

public class EmbeddedNumberField<T extends Element<T>, F extends ValueField<T>> extends AbstractField<NFE>
		implements ValueField<NFE> {
	private NumberField field;
	private F embeddingField;
	private T alpha;
	private MathMap<T, NFE> round;
	private boolean realEmbedding;

	EmbeddedNumberField(NumberField field, F embeddingField, T alpha, MathMap<T, NFE> round) {
		this.field = field;
		this.embeddingField = embeddingField;
		this.alpha = alpha;
		this.realEmbedding = alpha instanceof Real;
		this.round = round;
	}

	@Override
	public Exactness exactness() {
		return Exactness.EXACT;
	}

	@Override
	public boolean isComplete() {
		return false;
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

	public F embeddingField() {
		return embeddingField;
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
	public FactorizationResult<Polynomial<NFE>, NFE> factorization(UnivariatePolynomial<NFE> t) {
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

	public T embedding(NFE t) {
		Polynomial<T> asPolynomial = embeddingField.getUnivariatePolynomialRing().getEmbedding(t.asPolynomial(),
				new MathMap<>() {
					@Override
					public T evaluate(Fraction q) {
						return embeddingField.getFraction(q);
					}
				});
		return embeddingField.getUnivariatePolynomialRing().evaluate(asPolynomial, Collections.singletonList(alpha));
	}

	public MathMap<NFE, T> embeddingMap() {
		return new MathMap<>() {
			@Override
			public T evaluate(NFE t) {
				return embedding(t);
			}
		};
	}

	public NFE roundToInteger(T t) {
		return round.evaluate(t);
	}

	@Override
	public Real value(NFE t) {
		return embeddingField.value(embedding(t));
	}

	@Override
	public Reals getReals() {
		return embeddingField.getReals();
	}

	@Override
	public String toString() {
		return field.toString() + " in " + embeddingField.toString() + ", alpha=" + alpha;
	}
}
