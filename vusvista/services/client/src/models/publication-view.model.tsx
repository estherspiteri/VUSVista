export interface IPublicationPreview {
  publicationId: number;
  pmid?: string | null;
  title?: string | null;
  date?: Date | null;
  authors?: string[] | null;
  journal?: string | null;
  abstract?: string | null;
  isSupplementaryMaterialMatch: boolean;
  doi?: string | null;
  link?: string | null;
  isAddedManually: boolean;
}
