export interface IPublicationUpdates {
  title: string;
  doi?: string | null;
  link?: string | null;
}

export interface IVariantPublicationUpdates {
  lastEval: string;
  publicationUpdates: IPublicationUpdates[];
}
