package varieties.projective;

import java.util.ArrayList;
import java.util.List;

import fields.helper.CoordinateRing;
import fields.helper.CoordinateRing.CoordinateRingElement;
import fields.helper.FieldOfFractions.Fraction;
import fields.helper.LocalizedCoordinateRing;
import fields.helper.LocalizedCoordinateRing.LocalizedElement;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.local.Value;
import fields.local.ValueGroup;
import fields.polynomials.PolynomialIdeal;
import varieties.Morphism;
import varieties.affine.AffineMorphism;
import varieties.affine.AffinePoint;
import varieties.affine.AffineScheme;

public class ProjectiveMorphism<T extends Element<T>> implements Morphism<T, ProjectivePoint<T>, ProjectivePoint<T>> {
	private Field<T> field;
	private GenericProjectiveScheme<T> domain;
	private GenericProjectiveScheme<T> range;
	private List<Polynomial<T>> asPolynomials;
	private PolynomialRing<T> polynomials;

	public ProjectiveMorphism(GenericProjectiveScheme<T> domain, GenericProjectiveScheme<T> range,
			List<Polynomial<T>> asPolynomials) {
		if (!domain.getField().equals(range.getField())) {
			throw new ArithmeticException("Fields mismatched!");
		}
		this.field = domain.getField();
		this.domain = domain;
		this.range = range;
		this.asPolynomials = asPolynomials;
		int degree = -1;
		this.polynomials = domain.homogenousPolynomialRing();
		for (Polynomial<T> polynomial : asPolynomials) {
			if (polynomial.numberOfVariables() != polynomials.numberOfVariables()) {
				throw new ArithmeticException("domain dimensions mismatched!");
			}
			if (degree == -1) {
				degree = polynomial.degree();
			}
			if (!polynomials.isHomogeneous(polynomial)) {
				throw new ArithmeticException("polynomials not homogenous");
			}
			if (degree != polynomial.degree()) {
				throw new ArithmeticException("polynomials not of equal degree");
			}
		}
		if (asPolynomials.size() != range.homogenousPolynomialRing().numberOfVariables()) {
			throw new ArithmeticException("range dimensions mismatched!");
		}
	}
	
	public List<Polynomial<T>> asPolynomials() {
		return asPolynomials;
	}

	private class UpToUniformizerResult {
		private List<LocalizedElement<T>> upToUniformizer;
		private List<Value> values;
		private Value minimumValue;
	}

	private UpToUniformizerResult upToUniformizer(int coverIndex, PolynomialRing<T> affineRing,
			LocalizedCoordinateRing<T> localizedRing) {
		int[] map = new int[domain.homogenousPolynomialRing().numberOfVariables()];
		for (int i = 0; i < map.length; i++) {
			if (i == coverIndex) {
				map[i] = -1;
				continue;
			}
			map[i] = i > coverIndex ? i - 1 : i;
		}
		UpToUniformizerResult result = new UpToUniformizerResult();
		result.upToUniformizer = new ArrayList<>();
		result.values = new ArrayList<>();
		result.minimumValue = Value.INFINITY;
		for (Polynomial<T> polynomial : asPolynomials) {
			Polynomial<T> affinePolynomial = affineRing
					.getEmbedding(domain.homogenousPolynomialRing().dehomogenize(polynomial, coverIndex), map);
			LocalizedElement<T> inLocalization = localizedRing.getEmbedding(affinePolynomial);
			Value value = localizedRing.valuation(inLocalization);
			if (result.minimumValue.compareTo(value) > 0) {
				result.minimumValue = value;
			}
			result.values.add(value);
			result.upToUniformizer.add(localizedRing.upToUniformizerPower(inLocalization));
		}
		return result;
	}

	@Override
	public ProjectivePoint<T> evaluate(ProjectivePoint<T> t) {
		int coverIndex = domain.affineCoverIndex(t).get(0);
		AffineScheme<T> affineDomain = domain.getAffineCover().getCover().get(coverIndex);
		PolynomialRing<T> affineRing = affineDomain.getCoordinateRing().getPolynomialRing();
		AffinePoint<T> affinePreimage = domain.asAffinePoint(t, coverIndex);
		LocalizedCoordinateRing<T> localizedRing = affineDomain.localizedCoordinateRing(affinePreimage);
		UpToUniformizerResult result = upToUniformizer(coverIndex, affineRing, localizedRing);
		List<T> coords = new ArrayList<>();
		for (int i = 0; i < result.values.size(); i++) {
			if (result.values.get(i).equals(result.minimumValue)) {
				coords.add(localizedRing.reduce(result.upToUniformizer.get(i)));
			} else {
				coords.add(field.zero());
			}
		}
		return new ProjectivePoint<>(field, coords);
	}

	@Override
	public GenericProjectiveScheme<T> getDomain() {
		return domain;
	}

	@Override
	public GenericProjectiveScheme<T> getRange() {
		return range;
	}

	@Override
	public RestrictionResult<T> restrict(ProjectivePoint<T> preimage) {
		ProjectivePoint<T> image = evaluate(preimage);
		int domainCoverIndex = domain.affineCoverIndex(preimage).get(0);
		int rangeCoverIndex = range.affineCoverIndex(image).get(0);
		AffineScheme<T> affineDomain = domain.getAffineCover().getCover().get(domainCoverIndex);
		AffineScheme<T> affineRange = range.getAffineCover().getCover().get(rangeCoverIndex);
		PolynomialRing<T> affineRing = affineDomain.getCoordinateRing().getPolynomialRing();
		AffinePoint<T> affinePreimage = domain.asAffinePoint(preimage, domainCoverIndex);
		LocalizedCoordinateRing<T> localizedRing = affineDomain.localizedCoordinateRing(affinePreimage);
		UpToUniformizerResult upToUniformizer = upToUniformizer(domainCoverIndex, affineRing, localizedRing);
		if (!upToUniformizer.minimumValue.equals(upToUniformizer.values.get(rangeCoverIndex))) {
			throw new ArithmeticException("Range Cover Index wrong!");
		}
		Fraction<Polynomial<T>> rangeCoverFraction = upToUniformizer.upToUniformizer.get(rangeCoverIndex)
				.asPolynomialFraction();
		rangeCoverFraction.canonicalize();
		LocalizedElement<T> localDenominator = localizedRing.getEmbedding(rangeCoverFraction.getDenominator());

		Polynomial<T> denominator = rangeCoverFraction.getNumerator();
		List<Polynomial<T>> numerators = new ArrayList<>();
		List<Polynomial<T>> denominators = new ArrayList<>();
		ValueGroup g = ValueGroup.g();
		for (int i = 0; i < upToUniformizer.values.size(); i++) {
			if (i == rangeCoverIndex) {
				continue;
			}
			if (upToUniformizer.values.get(i).isInfinite()) {
				numerators.add(affineRing.zero());
				continue;
			}
			LocalizedElement<T> reduced = localizedRing.multiply(localDenominator,
					upToUniformizer.upToUniformizer.get(i), localizedRing.power(localizedRing.uniformizer(),
							g.subtract(upToUniformizer.values.get(i), upToUniformizer.minimumValue).value()));
			Fraction<Polynomial<T>> asFraction = reduced.asPolynomialFraction();
			asFraction.canonicalize();
			denominator = affineRing.lcm(asFraction.getDenominator(), denominator);
			numerators.add(asFraction.getNumerator());
			denominators.add(asFraction.getDenominator());
		}
		Polynomial<T> inverse = polynomials.subtract(polynomials.multiply(polynomials.getEmbedding(denominator),
				polynomials.getVar(polynomials.numberOfVariables())), polynomials.one());
		List<Polynomial<T>> inverseList = new ArrayList<>();
		inverseList.add(inverse);
		for (Polynomial<T> generator : affineDomain.getCoordinateRing().getIdeal().generators()) {
			inverseList.add(polynomials.getEmbedding(generator));
		}
		PolynomialIdeal<T> inverseIdeal = polynomials.getIdeal(inverseList);
		CoordinateRing<T> domainSubsetRing = new CoordinateRing<>(polynomials, inverseIdeal);
		AffineScheme<T> domainSubset = new AffineScheme<>(field, domainSubsetRing);
		List<CoordinateRingElement<T>> embeddingPolynomials = new ArrayList<>();
		for (int i = 0; i < affineRing.numberOfVariables(); i++) {
			embeddingPolynomials.add(domainSubsetRing.getEmbedding(polynomials.getVar(i + 1)));
		}
		AffineMorphism<T> embedding = new AffineMorphism<>(domainSubset, affineDomain, embeddingPolynomials);
		List<CoordinateRingElement<T>> asPolynomials = new ArrayList<>();
		for (int i = 0; i < numerators.size(); i++) {
			Polynomial<T> extend = affineRing.divideChecked(denominator, denominators.get(i));
			Polynomial<T> numerator = affineRing.multiply(numerators.get(i), extend);
			Polynomial<T> polynomial = polynomials.multiply(polynomials.getEmbedding(numerator),
					polynomials.getVar(polynomials.numberOfVariables()));
			asPolynomials.add(domainSubsetRing.getEmbedding(polynomial));
		}
		AffineMorphism<T> restricted = new AffineMorphism<>(domainSubset, affineDomain, asPolynomials);
		return new RestrictionResult<>(domainCoverIndex, rangeCoverIndex, embedding, affineRange.identityMorphism(),
				restricted);
	}

}
