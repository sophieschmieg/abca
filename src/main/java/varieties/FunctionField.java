package varieties;

import java.math.BigInteger;
import java.util.Iterator;

import fields.exceptions.InfinityException;
import fields.helper.AbstractField;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.polynomials.CoordinateRing;
import fields.polynomials.CoordinateRing.CoordinateRingElement;
import fields.polynomials.LocalizedCoordinateRing.LocalizedElement;
import varieties.affine.AffineScheme;
import varieties.curves.ProjectiveLine;
import varieties.projective.ProjectiveMorphism;
import varieties.projective.ProjectiveScheme;

public class FunctionField<T extends Element<T>> extends AbstractField<RationalFunction<T>> {
	private Field<T> field;
	private ProjectiveScheme<T> domain;
	private ProjectiveLine<T> range;
	private PolynomialRing<T> projectivePolynomialRing;
	private AffineScheme<T> affineSlice;
	private CoordinateRing<T> affineCoordinateRing;
	private PolynomialRing<T> affinePolynomialRing;
	private int affineCoverIndex;

	
	public FunctionField(ProjectiveScheme<T> domain) {
		this.domain = domain;
		this.field = domain.getField();
		this.projectivePolynomialRing = domain.asGenericProjectiveScheme().homogenousPolynomialRing();
		this.range = new ProjectiveLine<>(field);
		this.affineCoverIndex = domain.asGenericProjectiveScheme().homogenousPolynomialRing().numberOfVariables() - 1;
		this.affineSlice = domain.getAffineCover().getCover().get(affineCoverIndex);
		this.affineCoordinateRing = affineSlice.getCoordinateRing();
		this.affinePolynomialRing = affineCoordinateRing.getPolynomialRing();
	}

	@Override
	public Exactness exactness() {
		return field.exactness();
	}

	@Override
	public boolean isFinite() {
		return false;
	}

	public Field<T> getField() {
		return field;
	}

	public ProjectiveScheme<T> getDomain() {
		return domain;
	}

	public RationalFunction<T> getEmbedding(Polynomial<T> t) {
		if (t.equals(projectivePolynomialRing.zero())) {
			return getFunction(t, projectivePolynomialRing.one());
		}
		return getFunction(t, projectivePolynomialRing.getVarPower(affineCoverIndex + 1, t.degree()));
	}

	public RationalFunction<T> getEmbedding(T t) {
		return getEmbedding(projectivePolynomialRing.getEmbedding(t));
	}

	public RationalFunction<T> getFunction(Polynomial<T> numerator, Polynomial<T> denominator) {
		return RationalFunction.fromFraction(domain, range, numerator, denominator);
	}

	public RationalFunction<T> getFunction(CoordinateRingElement<T> numerator, CoordinateRingElement<T> denominator) {
		Polynomial<T> homogenousNumerator = domain.asGenericProjectiveScheme().homogenize(numerator.getElement(),
				affineCoverIndex);
		Polynomial<T> homogenousDenominator = domain.asGenericProjectiveScheme().homogenize(denominator.getElement(),
				affineCoverIndex);
		int degree = Math.max(homogenousNumerator.degree(), homogenousDenominator.degree());
		homogenousNumerator = projectivePolynomialRing.multiply(homogenousNumerator,
				projectivePolynomialRing.getVarPower(affineCoverIndex + 1, degree - homogenousNumerator.degree()));
		homogenousDenominator = projectivePolynomialRing.multiply(homogenousDenominator,
				projectivePolynomialRing.getVarPower(affineCoverIndex + 1, degree - homogenousDenominator.degree()));
		return getFunction(homogenousNumerator, homogenousDenominator);
	}
	
	public RationalFunction<T> getFunction(LocalizedElement<T> localizedElement) {
		return getFunction(localizedElement.getNumerator(), localizedElement.getDenominator());
	}
	
	public ProjectiveMorphism<T> asProjectiveMorphism(RationalFunction<T> t) {
		return t.asProjectiveMorphism();
	}
	
	public RationalFunction<T> fromProjectiveMorphism(ProjectiveMorphism<T> t) {
		if (t.asPolynomials().size() != 2) {
			throw new ArithmeticException("Not a projective morphism to the projective line!");
		}
		return getFunction(t.asPolynomials().get(0), t.asPolynomials().get(1));
	}
	
	public int affineCoverIndex() {
		return affineCoverIndex;
	}
	
	public AffineScheme<T> getAffineVariety() {
		return affineSlice;
	}
	
	public CoordinateRing<T> affineCoordinateRing() {
		return affineCoordinateRing;
	}

	@Override
	public RationalFunction<T> zero() {
		return getEmbedding(field.zero());
	}

	@Override
	public RationalFunction<T> one() {
		return getEmbedding(field.one());
	}

	@Override
	public BigInteger characteristic() {
		return this.field.characteristic();
	}
	private CoordinateRingElement<T> getAffine(Polynomial<T> t) {
		return affineCoordinateRing.getEmbedding(affinePolynomialRing.getEmbedding(projectivePolynomialRing.dehomogenize(t, affineCoverIndex+1)));
	}
		
	public CoordinateRingElement<T> getIntegerNumerator(RationalFunction<T> t) {
		return getAffine(t.getNumerator());
	}

	public CoordinateRingElement<T> getIntegerDenominator(RationalFunction<T> t) {
		return getAffine(t.getDenominator());
	}

	@Override
	public RationalFunction<T> add(RationalFunction<T> t1, RationalFunction<T> t2) {
		CoordinateRingElement<T> numerator = affineCoordinateRing.add(affineCoordinateRing.multiply(getIntegerNumerator(t1), getIntegerDenominator(t2)),
				affineCoordinateRing.multiply(getIntegerNumerator(t2), getIntegerDenominator(t1)));
		CoordinateRingElement<T> denominator = affineCoordinateRing.multiply(getIntegerDenominator(t1), getIntegerDenominator(t2));
		return getFunction(numerator, denominator);
	}

	@Override
	public RationalFunction<T> negative(RationalFunction<T> t) {
		return getFunction(projectivePolynomialRing.negative(t.getNumerator()), t.getDenominator());
	}

	@Override
	public RationalFunction<T> multiply(RationalFunction<T> t1, RationalFunction<T> t2) {
		CoordinateRingElement<T> numerator = affineCoordinateRing.multiply(getIntegerNumerator(t1), getIntegerNumerator(t2));
		CoordinateRingElement<T> denominator = affineCoordinateRing.multiply(getIntegerDenominator(t1), getIntegerDenominator(t2));
		return getFunction(numerator, denominator);
	}

	@Override
	public RationalFunction<T> inverse(RationalFunction<T> t) {
		return getFunction(t.getDenominator(), t.getNumerator());
	}

	@Override
	public RationalFunction<T> getRandomElement() {
		return null;
	}

	@Override
	public BigInteger getNumberOfElements() {
		return BigInteger.valueOf(-1);
	}

	@Override
	public Iterator<RationalFunction<T>> iterator() {
		throw new InfinityException();
	}
}
