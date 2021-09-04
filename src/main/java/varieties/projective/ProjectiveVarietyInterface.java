package varieties.projective;

import fields.interfaces.Element;
import varieties.ProjectivePoint;
import varieties.ProjectiveVariety;
import varieties.Variety;

public interface ProjectiveVarietyInterface<T extends Element<T>> extends Variety<T, ProjectivePoint<T>> {
	ProjectiveVariety<T> asProjectiveVariety();
}
