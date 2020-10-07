package varieties.curves;

import java.util.List;

import fields.helper.CoordinateRing;
import fields.helper.CoordinateRing.CoordinateRingElement;
import fields.interfaces.Element;
import fields.interfaces.Polynomial;
import varieties.ProjectivePoint;
import varieties.RationalFunction;
import varieties.Variety;
import varieties.curves.DivisorGroup.Divisor;

public interface SmoothCurve<T extends Element<T>> extends Variety<T> {
	public int getEmbeddingDimension();
	public boolean isProjective();
	public List<Polynomial<T>> getCotangentSpace(ProjectivePoint<T> p);
	public List<Polynomial<T>> getTangentSpace(ProjectivePoint<T> p);
	public List<RationalFunction<T>> getRiemannRochSpace(Divisor<T> div);
	public boolean isPrincipal(Divisor<T> div);
	public CoordinateRing<T> getCoordinateRing();
	public boolean hasSimplify();
	public List<CoordinateRingElement<T>> simplify(RationalFunction<T> t);
}
