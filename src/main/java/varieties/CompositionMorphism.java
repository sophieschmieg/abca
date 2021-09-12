package varieties;

import fields.interfaces.Element;
import varieties.affine.AffineMorphism;
import varieties.affine.AffineScheme;
import varieties.affine.AffineScheme.IntersectionResult;

public class CompositionMorphism<T extends Element<T>, S extends Element<S>, U extends Element<U>, V extends Element<V>>
		implements Morphism<T, S, V> {
	private Morphism<T, S, U> firstMorphism;
	private Morphism<T, U, V> secondMorphism;

	public CompositionMorphism(Morphism<T, S, U> firstMorphism, Morphism<T, U, V> secondMorphism) {
		this.firstMorphism = firstMorphism;
		this.secondMorphism = secondMorphism;
	}

	@Override
	public V evaluate(S t) {
		return secondMorphism.evaluate(firstMorphism.evaluate(t));
	}

	@Override
	public Scheme<T, S> getDomain() {
		return firstMorphism.getDomain();
	}

	@Override
	public Scheme<T, V> getRange() {
		return secondMorphism.getRange();
	}

	@Override
	public RestrictionResult<T> restrict(S preimage) {
		U intermediate = firstMorphism.evaluate(preimage);
		RestrictionResult<T> firstRestriction = firstMorphism.restrict(preimage);
		RestrictionResult<T> secondRestriction = secondMorphism.restrict(intermediate);
		AffineScheme<T> firstRange = firstMorphism.getRange().getAffineCover().getCover().get(firstRestriction.getRangeCoverIndex());
		AffineScheme<T> secondDomain = secondMorphism.getDomain().getAffineCover().getCover().get(secondRestriction.getDomainCoverIndex());
		IntersectionResult<T> intermediateIntersection = AffineScheme.intersect(firstRange, secondDomain);
		AffineMorphism<T> preimageMorphism = firstRestriction.getRestrictedMorphism().preimage(intermediateIntersection.getFirstEmbedding());
		AffineMorphism<T> restriction = AffineMorphism.composition(AffineMorphism.composition(preimageMorphism, firstRestriction.getRestrictedMorphism()), secondRestriction.getRestrictedMorphism());
		return new RestrictionResult<>(firstRestriction.getDomainCoverIndex(), secondRestriction.getRangeCoverIndex(), preimageMorphism, secondRestriction.getRangeEmbedding(), restriction);
	}
}
