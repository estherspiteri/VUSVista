export interface IPublicationPreview {
  pmid: string;
  title: string;
  date: Date;
  authors: string[];
  journal: string;
  abstract: string;
  isSupplementaryMaterialMatch: boolean;
}

export interface IPublicationSearch {
  isLitvarIdFound: boolean;
  publications: IPublicationPreview[];
}
