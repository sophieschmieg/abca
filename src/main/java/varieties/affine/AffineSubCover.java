package varieties.affine;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import fields.interfaces.Element;
import varieties.affine.AffineMorphism.AffineFiberedProduct;

public class AffineSubCover<T extends Element<T>> {
	private AffineCover<T> coarse;
	private AffineCover<T> fine;
	private List<AffineMorphism<T>> refinement;
	private List<Integer> coarseCoverIndex;

	public AffineSubCover(AffineCover<T> coarse, List<List<AffineMorphism<T>>> refinement) {
		this.coarse = coarse;
		this.refinement = new ArrayList<>();
		this.coarseCoverIndex = new ArrayList<>();
		List<AffineScheme<T>> fineCover = new ArrayList<>();
		List<List<AffineMorphism<T>>> fineAtlas = new ArrayList<>();
		for (int i = 0; i < coarse.getCover().size(); i++) {
			for (int j = 0; j < refinement.get(i).size(); j++) {
				AffineMorphism<T> morphism = refinement.get(i).get(j);
				if (!morphism.getRange().equals(coarse.getCover().get(i))) {
					throw new ArithmeticException("Not a subcover!");
				}
				if (!morphism.isInjective()) {
					throw new ArithmeticException("Not an embedding!");
				}
				fineCover.add(morphism.getDomain());
				this.refinement.add(morphism);
				this.coarseCoverIndex.add(i);
				List<AffineMorphism<T>> intersections = new ArrayList<>();
				for (int k = 0; k < coarse.getCover().size(); k++) {
					AffineMorphism<T> coarseIntersection = coarse.getIntersectionEmbedding(i, k);
					AffineFiberedProduct<T> refinementCoarseIntersection = AffineMorphism
							.fiberedProduct(coarseIntersection, morphism);
					AffineMorphism<T> toCoarseIntersection = refinementCoarseIntersection.getProjection1();
					AffineMorphism<T> toOtherCoarse = AffineMorphism.composition(toCoarseIntersection,
							coarse.getIntersectionEmbedding(k, i));
					AffineMorphism<T> toRefinement = refinementCoarseIntersection.getProjection2();
					for (int l = 0; l < refinement.get(k).size(); l++) {
						AffineMorphism<T> otherMorphism = refinement.get(k).get(l);
						AffineFiberedProduct<T> otherIntersection = AffineMorphism.fiberedProduct(toOtherCoarse,
								otherMorphism);
						intersections.add(AffineMorphism.composition(otherIntersection.getProjection2(), toRefinement));
					}
				}
				fineAtlas.add(intersections);
			}
		}
		this.fine = new AffineCover<>(fineCover, fineAtlas);
	}

	private static <T extends Element<T>> List<List<AffineMorphism<T>>> trivialRefinement(AffineCover<T> cover) {
		List<List<AffineMorphism<T>>> result = new ArrayList<>();
		for (int i = 0; i < cover.getCover().size(); i++) {
			result.add(Collections.singletonList(cover.getCover().get(i).identityMorphism()));
		}
		return result;
	}

	public AffineSubCover(AffineCover<T> cover) {
		this(cover, trivialRefinement(cover));
	}

	public AffineCover<T> getCoarseCover() {
		return coarse;
	}

	public AffineCover<T> getFineCover() {
		return fine;
	}

	public AffineMorphism<T> getRefinement(int fineCoverIndex) {
		return refinement.get(fineCoverIndex);
	}

	public int getCoarseCoverIndex(int fineCoverIndex) {
		return coarseCoverIndex.get(fineCoverIndex);
	}
}
