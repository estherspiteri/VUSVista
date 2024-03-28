import { IPublicationPreview } from "../models/publication-view.model";

export function convertPubDates(
  publications: IPublicationPreview[]
): IPublicationPreview[] {
  let updatedPublications: IPublicationPreview[] = [];

  if (publications) {
    //format publication date
    publications.forEach((publication) => {
      let pub = publication;
      pub.date = new Date(publication.date);
      updatedPublications = updatedPublications.concat(pub);
    });
  }

  return updatedPublications;
}
