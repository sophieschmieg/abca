package fields.polynomials;

import java.io.IOException;
import java.math.BigInteger;
import java.util.Iterator;

import fields.exceptions.InfinityException;
import fields.helper.AbstractField;
import fields.helper.TranscendentalFieldExtension;
import fields.helper.TranscendentalFieldExtension.TExt;
import fields.interfaces.AlgebraicExtensionElement;
import fields.interfaces.DedekindRing;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.FieldExtension;
import fields.interfaces.GlobalField;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.UnivariatePolynomial;
import util.PeekableReader;

public class UnivariateTranscendentalFieldExtension<B extends Element<B>, E extends AlgebraicExtensionElement<B, E>, FE extends FieldExtension<B, E, FE>> extends AbstractField<TExt<E>>
		implements GlobalField<TExt<E>, Polynomial<E>, E> {
	private TranscendentalFieldExtension<E> field;
	private FE base;
	private FieldExtensionUnivariatePolynomialRing<B, E, FE> integers;

	public UnivariateTranscendentalFieldExtension( FE base) {
		this.base = base;
		this.field = new TranscendentalFieldExtension<>(base, new String[] { "X" });
	}

	public String toString() {
		return field.toString();
	}

	public TExt<E> parse(PeekableReader reader) throws IOException {
		return field.parse(reader);
	}

	public TExt<E> zero() {
		return field.zero();
	}

	public TExt<E> add(TExt<E> s1, TExt<E> s2) {
		return field.add(s1, s2);
	}

	public TExt<E> negative(TExt<E> s) {
		return field.negative(s);
	}

	public Exactness exactness() {
		return base.exactness();
	}

	public TExt<E> getRandomElement() {
		return field.getRandomElement();
	}

	public boolean isFinite() {
		return false;
	}

	public BigInteger getNumberOfElements() throws InfinityException {
		throw new InfinityException();
	}

	public Iterator<TExt<E>> iterator() {
		return field.iterator();
	}

	public TExt<E> one() {
		return field.one();
	}

	public BigInteger characteristic() {
		return base.characteristic();
	}

	public TExt<E> multiply(TExt<E> t1, TExt<E> t2) {
		return field.multiply(t1, t2);
	}

	public TExt<E> inverse(TExt<E> t) {
		return field.inverse(t);
	}

	public TExt<E> getEmbedding(E t) {
		return field.getEmbedding(t);
	}

	public MathMap<E, TExt<E>> getEmbeddingMap() {
		return field.getEmbeddingMap();
	}

	public TExt<E> getEmbedding(Polynomial<E> t) {
		return field.getEmbedding(t);
	}

	public TExt<E> getElement(Polynomial<E> numerator, Polynomial<E> denomiator) {
		return field.getElement(numerator, denomiator);
	}

	public <S extends Element<S>> TExt<E> getEmbedding(TExt<S> t, MathMap<S, E> map) {
		return field.getEmbedding(t, map);
	}

	public TExt<E> getVar() {
		return field.getVar(1);
	}

	public TExt<E> getVarPower(int power) {
		return field.getVarPower(1, power);
	}

	public Field<E> getBaseField() {
		return base;
	}

	public PolynomialRing<E> polynomialRing() {
		return field.polynomialRing();
	}

	public TExt<E> characteristicRoot(TExt<E> t, int power) {
		return field.characteristicRoot(t, power);
	}

	public boolean hasCharacteristicRoot(TExt<E> t, int power) {
		return field.hasCharacteristicRoot(t, power);
	}

	@Override
	public FactorizationResult<Polynomial<TExt<E>>, TExt<E>> factorization(UnivariatePolynomial<TExt<E>> t) {
		return field.factorization(t);
	}

	@Override
	public TExt<E> getInteger(Polynomial<E> t) {
		return getEmbedding(t);
	}

	@Override
	public FieldExtensionUnivariatePolynomialRing<B, E, FE> ringOfIntegers() {
		if (integers == null) {
			integers = new FieldExtensionUnivariatePolynomialRing<>(base);
		}
		return integers;
	}

	@Override
	public ExtensionOfGlobalField<TExt<E>, Polynomial<E>, E, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?> getGlobalFieldExtension(
			UnivariatePolynomial<TExt<E>> minimalPolynomial) {
		// TODO Auto-generated method stub
		return null;
	}
}
