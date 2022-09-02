package varieties;

import java.util.ArrayList;
import java.util.List;
import java.util.Optional;

import fields.helper.TranscendentalFieldExtension;
import fields.helper.TranscendentalFieldExtension.TExt;
import fields.interfaces.Element;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.local.Value;
import fields.polynomials.CoordinateRing;
import fields.polynomials.CoordinateRing.CoordinateRingElement;
import fields.polynomials.LocalizedCoordinateRing;
import fields.polynomials.LocalizedCoordinateRing.LocalizedElement;
import varieties.affine.AffineMorphism;
import varieties.affine.AffineMorphism.AffineFiberedProduct;
import varieties.affine.AffineMorphism.ExpansionResult;
import varieties.affine.AffinePoint;

public class GeneralRationalFunction<T extends Element<T>, S extends Element<S>, U extends Element<U>> {
	private Scheme<T, S> domain;
	private Scheme<T, U> range;
	private AffineMorphism<T> morphism;
	private int domainCoverIndex;
	private AffineMorphism<T> domainEmbedding;
	private int rangeCoverIndex;

	public GeneralRationalFunction(Scheme<T, S> domain, Scheme<T, U> range, AffineMorphism<T> morphism,
			int domainCoverIndex, AffineMorphism<T> domainEmbedding, int rangeCoverIndex) {
		this.domain = domain;
		this.range = range;
		this.domainCoverIndex = domainCoverIndex;
		this.domainEmbedding = domainEmbedding;
		this.morphism = morphism;
		this.rangeCoverIndex = rangeCoverIndex;
		if (!morphism.getDomain().equals(domainEmbedding.getDomain())) {
			throw new ArithmeticException("Domains not equal");
		}
		if (!domainEmbedding.getRange().equals(domain.getAffineCover().get(domainCoverIndex))) {
			throw new ArithmeticException("Domain embedding not equal");
		}
		if (!domainEmbedding.isOpenImmersion()) {
			throw new ArithmeticException("Domain not an open immersion!");
		}
	}

	public Scheme<T, S> getDomain() {
		return domain;
	}

	public Scheme<T, U> getRange() {
		return range;
	}

	public AffineMorphism<T> getMorphism() {
		return morphism;
	}

	public int getDomainCoverIndex() {
		return domainCoverIndex;
	}

	public AffineMorphism<T> getDomainEmbedding() {
		return domainEmbedding;
	}

	public int getRangeCoverIndex() {
		return rangeCoverIndex;
	}

	public TExt<T> inducedMap(CoordinateRingElement<T> t, int rangeCoverIndex, int domainCoverIndex) {
		CoordinateRing<T> domainCoordinateRing = domain.getAffineCover().get(domainCoverIndex).getCoordinateRing();
		PolynomialRing<T> domainPolynomialRing = domainCoordinateRing.getPolynomialRing();
		TranscendentalFieldExtension<T> domainTranscendental = new TranscendentalFieldExtension<>(domain.getField(),
				domainPolynomialRing);
		CoordinateRing<T> rangeCoordinateRing = range.getAffineCover().get(this.rangeCoverIndex).getCoordinateRing();
		if (this.domainCoverIndex == domainCoverIndex && this.rangeCoverIndex == rangeCoverIndex) {
			CoordinateRingElement<T> inDomainRestriction = morphism.inducedMap(t);
			return domainEmbedding.birationalInverseInducedMap(inDomainRestriction);
		}
		if (this.domainCoverIndex == domainCoverIndex) {
			AffineMorphism<T> rangeIntersection = range.getAffineCover().getIntersectionEmbedding(rangeCoverIndex,
					this.rangeCoverIndex);
			CoordinateRingElement<T> inIntersection = rangeIntersection.inducedMap(t);
			TExt<T> inRangeCover = range.getAffineCover()
					.getIntersectionEmbedding(this.rangeCoverIndex, rangeCoverIndex).birationalInverse()
					.inducedMap(inIntersection, 0, 0);
			TExt<T> numerator = inducedMap(rangeCoordinateRing.getEmbedding(inRangeCover.getNumerator()),
					this.rangeCoverIndex, domainCoverIndex);
			TExt<T> denominator = inducedMap(rangeCoordinateRing.getEmbedding(inRangeCover.getDenominator()),
					this.rangeCoverIndex, domainCoverIndex);
			return domainTranscendental.divide(numerator, denominator);
		}
		TExt<T> inDomainCover = inducedMap(t, rangeCoverIndex, this.domainCoverIndex);
		AffineMorphism<T> domainIntersection = domain.getAffineCover().getIntersectionEmbedding(this.domainCoverIndex,
				domainCoverIndex);
		Polynomial<T> inIntersectionNumerator = domainIntersection.inducedMap(inDomainCover.getNumerator());
		Polynomial<T> inIntersectionDenominator = domainIntersection.inducedMap(inDomainCover.getDenominator());
		GeneralRationalFunction<T, AffinePoint<T>, AffinePoint<T>> domainInverse = domain.getAffineCover()
				.getIntersectionEmbedding(domainCoverIndex, this.domainCoverIndex).birationalInverse();
		TExt<T> numerator = domainInverse.inducedMap(
				domainIntersection.getDomain().getCoordinateRing().getEmbedding(inIntersectionNumerator), 0, 0);
		TExt<T> denominator = domainInverse.inducedMap(
				domainIntersection.getDomain().getCoordinateRing().getEmbedding(inIntersectionDenominator), 0, 0);
		return domainTranscendental.divide(numerator, denominator);
	}

	public Optional<U> evaluate(S point) {
		if (domain.affineCoverIndex(point).contains(domainCoverIndex)) {
			AffinePoint<T> affine = domain.asAffinePoint(point, domainCoverIndex);
			List<AffinePoint<T>> preimage = domainEmbedding.preimageList(affine);
			if (preimage.size() == 1) {
				return Optional.of(range.fromAffinePoint(preimage.get(0), rangeCoverIndex));
			}
		}
		int domainIndex = domain.affineCoverIndex(point).get(0);
		AffinePoint<T> affine = domain.asAffinePoint(point, domainIndex);
		LocalizedCoordinateRing<T> localized = null;
		if (domain.dimension() == 1 && !domain.getAffineCover().get(domainIndex).singularPoints().contains(affine)) {
			CoordinateRing<T> cr = domain.getAffineCover().get(domainIndex).getCoordinateRing();
			localized = new LocalizedCoordinateRing<T>(domain.getField(), cr,
					cr.getIdeal(affine.asIdeal(cr.getPolynomialRing())));
		}
		AffineMorphism<T> transition = domain.getAffineCover().getIntersectionEmbedding(domainCoverIndex, domainIndex);
		AffineFiberedProduct<T> restrictedDomain = AffineMorphism.fiberedProduct(domainEmbedding, transition);
		AffineMorphism<T> movedMorphism = AffineMorphism.composition(restrictedDomain.getSimplifiedProjection1(),
				morphism);
		AffineMorphism<T> movedEmbedding = AffineMorphism.composition(restrictedDomain.getSimplifiedProjection2(),
				domain.getAffineCover().getIntersectionEmbedding(domainIndex, domainCoverIndex));
		rangeCoverLoop: for (int i = 0; i < range.getAffineCover().size(); i++) {
			AffineFiberedProduct<T> fiber = AffineMorphism.fiberedProduct(movedMorphism,
					range.getAffineCover().getIntersectionEmbedding(rangeCoverIndex, i));
			AffineMorphism<T> currentMorphism = AffineMorphism.composition(fiber.getSimplifiedProjection2(),
					range.getAffineCover().getIntersectionEmbedding(i, rangeCoverIndex));
			AffineMorphism<T> currentDomainEmbedding = AffineMorphism.composition(fiber.getSimplifiedProjection1(),
					movedEmbedding);
			ExpansionResult<T> currentExpansion = currentMorphism.expand(currentDomainEmbedding);
			List<AffinePoint<T>> currentPreimage = currentExpansion.getExpandedOpenImmersion().preimageList(affine);
			if (currentPreimage.size() != 1) {
				if (localized != null) {
					List<T> result = new ArrayList<>();
					for (Polynomial<T> f : currentExpansion.getExpandedMorphism().getPolynomials()) {
						TExt<T> inverse = currentExpansion.getExpandedOpenImmersion().birationalInverseInducedMap(currentExpansion.getExpandedOpenImmersion().getDomain()
								.getCoordinateRing().getEmbedding(f));
						LocalizedElement<T> asLocalizedElement = localized.getEmbedding(inverse.getNumerator(),
								inverse.getDenominator());
						Value v = localized.valuation(asLocalizedElement);
						if (v.compareTo(Value.ZERO) < 0) {
							continue rangeCoverLoop;
						}
						result.add(localized.reduceInteger(asLocalizedElement));
					}
					return Optional.of(range.fromAffinePoint(new AffinePoint<>(range.getField(), result), i));
				}
				continue;
			}
			AffinePoint<T> evaluated = currentExpansion.getExpandedMorphism().evaluate(currentPreimage.get(0));
			return Optional.of(range.fromAffinePoint(evaluated, i));
		}
		return Optional.empty();
	}
}
