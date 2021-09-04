package varieties;

import fields.interfaces.Element;
import fields.interfaces.MathMap;
import varieties.affine.AffineMorphism;

public interface Morphism<T extends Element<T>, S extends Element<S>, U extends Element<U>> extends MathMap<S, U> {
	Variety<T, S> getDomain();

	Variety<T, U> getRange();

	public class RestrictionResult<T extends Element<T>> {
		private int domainCoverIndex;
		private int rangeCoverIndex;
		private AffineMorphism<T> domainEmbedding;
		private AffineMorphism<T> rangeEmbedding;
		private AffineMorphism<T> restrictedMorphism;

		public RestrictionResult(int domainCoverIndex, int rangeCoverIndex, AffineMorphism<T> domainEmbedding,
				AffineMorphism<T> rangeEmbedding, AffineMorphism<T> restrictedMorphism) {
			this.domainCoverIndex = domainCoverIndex;
			this.rangeCoverIndex = rangeCoverIndex;
			this.domainEmbedding = domainEmbedding;
			this.rangeEmbedding = rangeEmbedding;
			this.restrictedMorphism = restrictedMorphism;
		}

		public int getDomainCoverIndex() {
			return domainCoverIndex;
		}

		public int getRangeCoverIndex() {
			return rangeCoverIndex;
		}

		public AffineMorphism<T> getDomainEmbedding() {
			return domainEmbedding;
		}

		public AffineMorphism<T> getRangeEmbedding() {
			return rangeEmbedding;
		}

		public AffineMorphism<T> getRestrictedMorphism() {
			return restrictedMorphism;
		}
	}

	RestrictionResult<T> restrict(S preimage);
}
