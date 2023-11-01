import { IPublicationSearch } from "../../models/publication-search/publication-search.model";

export interface IGetPublicationsRequest {
  rsid: string;
}

export interface IGetPublicationsResponse {
  isSuccess: boolean;
  publicationSearch?: IPublicationSearch;
}
