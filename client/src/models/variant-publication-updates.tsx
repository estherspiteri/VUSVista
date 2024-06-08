export interface IPublicationUpdates {
  title: string;
  doi?: string | null;
  link?: string | null;
  isManuallyAdded: boolean;
}

export interface IVariantPublicationUpdates {
  lastEval: string;
  publicationUpdates: IPublicationUpdates[];
}
