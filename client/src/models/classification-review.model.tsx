export interface IClassificationReviewPublications {
  title: string;
  doi?: string;
  link?: string;
}

export interface IClassificationReview {
  classification: string;
  reason: string;
  dateAdded: Date;
  scientificMemberName: string;
  scientificMemberEmail: string;
  acmgRules: string[];
  publications: IClassificationReviewPublications[];
}
