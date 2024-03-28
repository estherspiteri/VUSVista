export interface IPublicationPreview {
  publicationId: number;
  pmid?: string | null;
  title?: string | null;
  date?: Date | null; //TODO: fix to correct camelCase - change mapping from backend
  authors?: string[] | null;
  journal?: string | null;
  abstract?: string | null;
  isSupplementaryMaterialMatch: boolean;
  doi?: string | null;
  link?: string | null;
}

