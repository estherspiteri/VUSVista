import { IClassificationReview } from "../models/classification-review.model";
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

export function convertReviewDates(
  reviews: IClassificationReview[]
): IClassificationReview[] {
  let updatedReviews: IClassificationReview[] = [];

  if (reviews) {
    //format review date
    reviews.forEach((review) => {
      let r = review;
      r.dateAdded = new Date(review.dateAdded);
      updatedReviews = updatedReviews.concat(r);
    });
  }

  return updatedReviews;
}

export function getMonthString(month: number) {
  const months = [
    "Jan",
    "Feb",
    "Mar",
    "Apr",
    "May",
    "Jun",
    "Jul",
    "Aug",
    "Sep",
    "Oct",
    "Nov",
    "Dec",
  ];

  return months[month];
}
